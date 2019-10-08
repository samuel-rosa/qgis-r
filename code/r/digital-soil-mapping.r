##Modelling=group
##Spatial environmental modelling=name
##Observations=vector
##Target=field Observations
##Weights=field Observations
##Validation=field Observations
##Covariates=multiple raster
##Learner=selection Classification and Regression Tree (C & R);Linear Discriminant Analysis (C);Linear Regression (R);Linear Regression with Stepwise Selection (R);Penalized Multinomial Regression (C);Neural Network (C & R);Random Forest (C & R);Support Vector Machines w/ Radial Basis Function Kernel (C & R)
##Predictions=output raster
##Uncertainty=output raster
##Metadata=output table

# Load necessary packages ----
library(sp)
library(raster)
library(caret)
library(snow)

# Check integrity of data type ----
if (is.factor(Observations[[Validation]])) {
  Observations[[Validation]] <- as.integer(levels(Observations[[Validation]]))[Observations[[Validation]]]
}
if (is.factor(Observations[[Weights]])) {
  Observations[[Weights]] <- as.numeric(levels(Observations[[Weights]]))[Observations[[Weights]]]
}

# Remove observations with NAs ----
na_idx <- complete.cases(Observations@data)
Observations <- Observations[na_idx, ]

# Identify validation observations (if any) ----
if (any(Observations[[Validation]]) == 1) {
  validate <- TRUE
  idx <- which(Observations[[Validation]] == 1)
  val_data <- Observations[idx, ]@data
  Observations <- Observations[-idx, -which(colnames(Observations@data) == Validation)]
  Observations[[Target]] <- as.factor(as.character(Observations[[Target]]))
  n_val <- length(val_data[[Target]])
} else {
  validate <- FALSE
}

# Identify learner and set arguments, including the type of spatial predicions ----
model <- c("rpart", "lda", "lm", "lmStepAIC", "multinom", "nnet", "rf", "svmRadial")
Learner <- Learner + 1
model <- model[Learner]
if (is.numeric(Observations[[Target]])) {
  index <- 1
  type <- "raw"
} else {
  index <- 1:nlevels(Observations[[Target]])
  type <- "prob"
}
if (model == "svmRadial" & type == "prob") {
  prob.model <- TRUE
} else {
  prob.model <- FALSE
}

# Funtion to compute the Shannon entropy ----
entropy <-
  function (x) {
    - sum(x * log(x, base = length(x)), na.rm = TRUE)
  }

# Function to compute the confusion index ----
confusion <-
  function (x) {
    1 - diff(sort(x, decreasing = TRUE)[2:1])
  }

# Look for factor covariates ----
# I am not sure if factor covariates must be specified as such.
covar_cols <- which(!colnames(Observations@data) %in% c(Target, Weights, Validation))
# if (any(sapply(Observations@data[, covar_cols], is.factor))) {
#   is_char <- which(sapply(Observations@data[, covar_cols], is.factor))
#   is_char <- match(names(is_char), sapply(Covariates, names))
#   for (i in 1:length(is_char)) {
#     if (!is.factor(Covariates[[is_char[i]]])) {
#       Covariates[[is_char[i]]] <- as.factor(Covariates[[is_char[i]]]) 
#     }
#   }
# }

# Calibrate statistical learner ----
# We must pass a formula to avoid the errors reported in https://stackoverflow.com/a/25272143/3365410 
form <- formula(paste(Target, " ~ ", paste(colnames(Observations@data[, covar_cols]), collapse = " + ")))
# rpart does not deal with parameters that are supposed to be passed to other models because they are checked
# against a list of valid function arguments.
if (model == "rpart") {
  learner_fit <- train(
    form = form, data = Observations@data, weights = Observations[[Weights]], method = model, tuneLength = 1,
    na.action = na.omit, trControl = trainControl(method = "LOOCV")
  )
} else {
  learner_fit <- train(
    form = form, data = Observations@data, weights = Observations[[Weights]], method = model, tuneLength = 1,
    na.action = na.omit, trControl = trainControl(method = "LOOCV"), 
    importance = TRUE, #rf
    prob.model = prob.model # svmRadial
  )
}

# Perform validation if validation data is available ----
if (validate) {
  pred <- predict(learner_fit, val_data)
  if (type == "raw") {
    # nothing defined yet
  } else {
    error <- confusionMatrix(data = pred, reference = val_data[[Target]])
  }
}

# Make spatial predictions ----
beginCluster()
prediction <- 
  raster::clusterR(
    raster::brick(Covariates), 
    raster::predict, 
    args = list(model = learner_fit, type = type, index = index)
  )
raster::endCluster()

# Compute predictions and prediction uncertainty ----
if (type == "prob") {
  Predictions <- as.factor(calc(x = prediction, fun = nnet::which.is.max))
  rat <- levels(Predictions)[[1]]
  rat$class <- levels(Observations[[Target]])[rat$ID]
  levels(Predictions) <- rat
  Uncertainty <-
    brick(
      calc(x = prediction, fun = max),
      calc(x = prediction, fun = entropy),
      calc(x = prediction, fun = confusion)
    )
  Metadata <- 
    rbind(
      c("Predictions", 
        paste(
          "Predicted class (", paste(apply(rat, 1, paste, collapse = "="), collapse = "; ", sep = ""), "); ",
          "Observations = ", nrow(learner_fit$trainingData), sep = "")),
      c("Uncertainty", 
        "Band 1 = Theoretical purity (0-1); Band 2 = Shannon entropy (0-1); Band 3 = Confusion index (0-1)"),
      c("Statistical learner", 
        paste(
          learner_fit$method[1], " = ", learner_fit$modelInfo$label[1], " (", learner_fit$modelType[1], ")",
          sep = "")),
      c("Cross-validation", 
        paste("Overall accuracy = ", round(learner_fit$results$Accuracy[nrow(learner_fit$results)], 4), "; ",
              "Overall kappa = ", round(learner_fit$results$Kappa[nrow(learner_fit$results)], 4), 
              sep = "")),
      c("Covariate importance",
        paste(rownames(varImp(learner_fit)[[1]]), collapse = "; "))
    )
  if (validate) {
    Metadata <- 
      rbind(
        Metadata,
        c("Validation", 
          paste("Overall accuracy = ", round(error$overall[["Accuracy"]], 4), "; ",
                "Overall kappa = ", round(error$overall[["Kappa"]], 4), "; ",
                "Observations = ", n_val,
                sep = ""))
      )
  }
  colnames(Metadata) <- c("Item", "Description")
} else { # Quantitative variable
  
  # Use predictions of nonlinear model as covariate in a linear model
  if (!model %in% c("lm", "lmStepAIC")) {
    Observations@data$prediction <- as.vector(predict(learner_fit, Observations@data))
    form <- formula(paste(Target, " ~ prediction"))
    learner_fit2 <- caret::train(
      form = form, data = Observations@data, weights = Observations[[Weights]], method = "lm",
      na.action = na.omit, trControl = trainControl(method = "LOOCV"))
    names(prediction) <- "prediction"
    
    # Fit linear regression model using lm() to be able to compute (approximate) prediction intervals 
    lm_fit <- lm(formula = form, data = Observations@data)
    beginCluster()
    prediction <-
      raster::clusterR(
        prediction,
        raster::predict,
        args = list(fun = predict.lm, model = lm_fit, interval = 'prediction', index = 1:3)
      )
    raster::endCluster()
  }
  names(prediction) <- c("fit", "lwr", "upr")
  Predictions <- prediction[[1]]
  Uncertainty <- prediction[[2:3]]
  
  Metadata <-
    rbind(
      c("Predictions",
        paste("Predicted values (", Target, "); ", "Observations = ", nrow(learner_fit$trainingData),
              sep = "")),
      c("Uncertainty",
        # "Band 1 = Lower prediction limit (95%); Band 2 = Upper prediction limit (95%); Band 3 = Prediction error standard deviation"),
        "Band 1 = Lower prediction limit (95%); Band 2 = Upper prediction limit (95%)"),
      # c("Uncertainty", paste("Predicted values (", Target, ")", sep = "")), # temporary
      c("Statistical learner",
        paste(learner_fit$method[1], " = ", learner_fit$modelInfo$label[1], " (", learner_fit$modelType[1], ")",
              sep = "")),
      c("Cross-validation",
        paste(
          'ME = ', round(mean(apply(learner_fit2$pred[c('pred', 'obs')], 1, diff)), 4), "; ",
          'MAE = ', round(learner_fit2$results$MAE[nrow(learner_fit2$results)], 4), "; ",
          "RMSE = ", round(learner_fit2$results$RMSE[nrow(learner_fit2$results)], 4), "; ",
          "AVE = ", round(1 - sum(apply(learner_fit2$pred[c('pred', 'obs')], 1, diff)^2)/sum((mean(learner_fit2$pred$obs) - learner_fit2$pred$obs)^2), 4), sep = "")),
      c("Covariate importance",
        paste(rownames(varImp(learner_fit)[[1]])[order(varImp(learner_fit)[[1]], decreasing = TRUE)],
              collapse = "; "))
    )
  colnames(Metadata) <- c("Item", "Description")
}

# Output ----
Predictions
Uncertainty
Metadata
