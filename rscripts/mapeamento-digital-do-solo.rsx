##Modelling=group
##Mapeamento digital do solo=name
##Observations=vector
##Response=field Observations
##Weights=optional field Observations
##Validation=optional field Observations
##Covariates=multiple raster
##Model=selection Classification and Regression Tree (Cat & Num);Linear Discriminant Analysis (Cat);Linear Regression (Num);Linear Regression with Stepwise Selection (Num);Penalized Multinomial Regression (Cat) ;
##Predictions=output raster
##Uncertainty=output raster
##Metadata=output table

# Print R session info ----
>cat('\n\n\n')
>sessionInfo()
>cat('\n\n\n')

# Check loaded observations ----
>cat('\n\n\nThese are the observations that you loaded:\n\n\n')
>print(Observations)

# Check loaded covariates ----
if (is.null(Covariates)) {
  stop ("\n\n\nIt looks like you forgot to load the covariates!\n\n\n")
}
>cat('\n\n\nThese are the covariates that you loaded:\n\n\n')
>print(Covariates)
>cat('\n\n\n')

# Load necessary packages ----
library(caret)
library(snow)

# Drop geometry and set response data type ----
Observations <- sf::st_drop_geometry(Observations)
if (is.character(Observations[[Response]])) {
  Observations[[Response]] <- as.factor(Observations[[Response]])
}

# Weights ----
# Check if Weights is NULL
if (is.null(Weights)) {
  Weights <- 'weights'
  Observations$weights <- as.double(1)
}

# Validation ----
# Check if Validation is NULL
# Identify validation observations (if any)
if (is.null(Validation)) {
  validate <- FALSE
  Validation <- 'validation'
} else if (any(Observations[[Validation]]) == 1) {
  validate <- TRUE
  idx <- which(Observations[[Validation]] == 1)
  val_data <- Observations[idx, ]
  Observations <- Observations[-idx, -which(colnames(Observations) == Validation)]
  n_val <- length(val_data[[Response]])
} else {
  validate <- FALSE
}

# Check if all covariates were loaded ----
covar_cols <- which(!colnames(Observations) %in% c(Response, Weights, Validation))
covar_names <- colnames(Observations)[covar_cols]
if (!all(covar_names %in% sapply(Covariates, names))) {
  covar_out <- covar_names[which(!covar_names %in% sapply(Covariates, names))]
  covar_out <- paste(covar_out, collapse = ', ')
  stop (paste("\n\n\nIt looks like you forgot to load the following covariates:", covar_out, "\n\n\n"))
}

# Remove observations with NAs ----
na_idx <- complete.cases(Observations)
Observations <- Observations[na_idx, ]

# Identify model and set arguments, including the type of spatial predicions ----
model <- c("rpart", "lda", "lm", "lmStepAIC", "multinom")
Model <- Model + 1
model <- model[Model]
if (is.numeric(Observations[[Response]])) {
  index <- 1
  type <- "raw"
} else {
  index <- 1:nlevels(Observations[[Response]])
  type <- "prob"
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
covar_cols <- which(!colnames(Observations) %in% c(Response, Weights, Validation))

# Calibrate statistical model ----
# We must pass a formula to avoid the errors reported in https://stackoverflow.com/a/25272143/3365410 
form <- formula(paste(Response, " ~ ", paste(colnames(Observations[, covar_cols]), collapse = " + ")))
model_fit <- caret::train(
  form = form,
  data = Observations,
  weights = Observations[[Weights]], 
  method = model,
  na.action = na.omit,
  tuneLength = 1,
  trControl = caret::trainControl(method = "LOOCV")
)

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
>ifelse(model %in% c("lm", "lmStepAIC", "multinom"), summary(model_fit$finalModel), print(model_fit$finalModel))
>cat('\n\n\n')
