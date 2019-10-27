##Modelling=group
##Digital soil mapping=name
##Observations=vector
##Response=field Observations
##Weights=optional field Observations
##Validation=optional field Observations
##Covariates=multiple raster
##Model=selection Classification and Regression Tree (C & R);Linear Discriminant Analysis (C);Linear Regression (R);Linear Regression with Stepwise Selection (R);Penalized Multinomial Regression (C);Neural Network (C & R);Random Forest (C & R);Support Vector Machines w/ Radial Basis Function Kernel (C & R)
##Predictions=output raster
##Uncertainty=output raster
##Metadata=output table

# Load necessary packages ----
library(sp)
library(caret)
library(snow)

# Remove observations with NAs ----
na_idx <- complete.cases(Observations@data)
Observations <- Observations[na_idx, ]

# Weights ----
# Check if Weights is NULL
# Check integrity of data type ----
if (is.null(Weights)) {
  Weights <- 'weights'
  Observations@data$weights <- as.double(1)
} else {
  if (is.factor(Observations[[Weights]])) {
    Observations[[Weights]] <- as.numeric(levels(Observations[[Weights]]))[Observations[[Weights]]]
  }
}

# Validation ----
# Check if Validation is NULL
# Check integrity of data type
# Identify validation observations (if any)
if (is.null(Validation)) {
  validate <- FALSE
} else {
  if (is.factor(Observations[[Validation]])) {
    Observations[[Validation]] <- as.integer(levels(Observations[[Validation]]))[Observations[[Validation]]]
  }
  if (any(Observations[[Validation]]) == 1) {
    validate <- TRUE
    idx <- which(Observations[[Validation]] == 1)
    val_data <- Observations[idx, ]@data
    Observations <- Observations[-idx, -which(colnames(Observations@data) == Validation)]
    Observations[[Response]] <- as.factor(as.character(Observations[[Response]]))
    n_val <- length(val_data[[Response]])
  } else {
    validate <- FALSE
  }
}

# Identify model and set arguments, including the type of spatial predicions ----
model <- c("rpart", "lda", "lm", "lmStepAIC", "multinom", "nnet", "rf", "svmRadial")
Model <- Model + 1
model <- model[Model]
if (is.numeric(Observations[[Response]])) {
  index <- 1
  type <- "raw"
} else {
  index <- 1:nlevels(Observations[[Response]])
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
covar_cols <- which(!colnames(Observations@data) %in% c(Response, Weights, Validation))
# covar_cols <- which(!colnames(Observations@data) %in% c(Response, Weights))
# if (any(sapply(Observations@data[, covar_cols], is.factor))) {
#   is_char <- which(sapply(Observations@data[, covar_cols], is.factor))
#   is_char <- match(names(is_char), sapply(Covariates, names))
#   for (i in 1:length(is_char)) {
#     if (!is.factor(Covariates[[is_char[i]]])) {
#       Covariates[[is_char[i]]] <- as.factor(Covariates[[is_char[i]]]) 
#     }
#   }
# }

# Calibrate statistical model ----
# We must pass a formula to avoid the errors reported in https://stackoverflow.com/a/25272143/3365410 
form <- formula(paste(Response, " ~ ", paste(colnames(Observations@data[, covar_cols]), collapse = " + ")))
# rpart does not deal with parameters that are supposed to be passed to other models because they are checked
# against a list of valid function arguments.
if (model == "rpart") {
  model_fit <- caret::train(
    form = form, data = Observations@data, weights = Observations[[Weights]], method = model, tuneLength = 1,
    na.action = na.omit, trControl = trainControl(method = "LOOCV")
  )
} else {
  model_fit <- caret::train(
    form = form, data = Observations@data, weights = Observations[[Weights]], method = model, tuneLength = 1,
    na.action = na.omit, trControl = trainControl(method = "LOOCV"), 
    importance = TRUE, #rf
    prob.model = prob.model # svmRadial
  )
}

# Perform validation if validation data is available ----
if (validate) {
  pred <- predict(model_fit, val_data)
  if (type == "raw") {
    # nothing defined yet
  } else {
    error <- caret::confusionMatrix(data = pred, reference = val_data[[Response]])
  }
}

# Make spatial predictions ----
raster::beginCluster()
prediction <- 
  raster::clusterR(
    raster::brick(Covariates), 
    raster::predict, 
    args = list(model = model_fit, type = type, index = index)
  )
raster::endCluster()

# Compute predictions and prediction uncertainty ----
if (type == "prob") {
  Predictions <- as.factor(calc(x = prediction, fun = nnet::which.is.max))
  rat <- levels(Predictions)[[1]]
  rat$class <- levels(Observations[[Response]])[rat$ID]
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
          "Observations = ", nrow(model_fit$trainingData), sep = "")),
      c("Uncertainty", 
        "Band 1 = Theoretical purity (0-1); Band 2 = Shannon entropy (0-1); Band 3 = Confusion index (0-1)"),
      c("Statistical model", 
        paste(
          model_fit$method[1], " = ", model_fit$modelInfo$label[1], " (", model_fit$modelType[1], ")",
          sep = "")),
      c("Cross-validation (overall)", 
        paste("Accuracy = ", round(model_fit$results$Accuracy[nrow(model_fit$results)], 4), "; ",
              "Kappa = ", round(model_fit$results$Kappa[nrow(model_fit$results)], 4), 
              sep = "")),
      c("Covariate importance",
        paste(rownames(caret::varImp(model_fit)[[1]]), collapse = "; "))
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
    Observations@data$prediction <- as.vector(predict(model_fit, Observations@data))
    form <- formula(paste(Response, " ~ prediction"))
    model_fit2 <- caret::train(
      form = form, data = Observations@data, weights = Observations[[Weights]], method = "lm",
      na.action = na.omit, trControl = caret::trainControl(method = "LOOCV"))
    names(prediction) <- "prediction"
    
    # Fit linear regression model using lm() to be able to compute (approximate) prediction intervals 
    lm_fit <- lm(formula = form, data = Observations@data)
    raster::beginCluster()
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
        paste("Predicted values (", Response, "); ", "Observations = ", nrow(model_fit$trainingData),
              sep = "")),
      c("Uncertainty",
        # "Band 1 = Lower prediction limit (95%); Band 2 = Upper prediction limit (95%); Band 3 = Prediction error standard deviation"),
        "Band 1 = Lower prediction limit (95%); Band 2 = Upper prediction limit (95%)"),
      # c("Uncertainty", paste("Predicted values (", Response, ")", sep = "")), # temporary
      c("Statistical model",
        paste(model_fit$method[1], " = ", model_fit$modelInfo$label[1], " (", model_fit$modelType[1], ")",
              sep = "")),
      c("Cross-validation",
        paste(
          'ME = ', round(mean(apply(model_fit2$pred[c('pred', 'obs')], 1, diff)), 4), "; ",
          'MAE = ', round(model_fit2$results$MAE[nrow(model_fit2$results)], 4), "; ",
          "RMSE = ", round(model_fit2$results$RMSE[nrow(model_fit2$results)], 4), "; ",
          "AVE = ", round(1 - sum(apply(model_fit2$pred[c('pred', 'obs')], 1, diff)^2)/sum((mean(model_fit2$pred$obs) - model_fit2$pred$obs)^2), 4), sep = "")),
      c("Covariate importance",
        paste(rownames(caret::varImp(model_fit)[[1]])[order(caret::varImp(model_fit)[[1]], decreasing = TRUE)],
              collapse = "; "))
    )
  colnames(Metadata) <- c("Item", "Description")
}

# Output ----
Predictions
Uncertainty
Metadata
