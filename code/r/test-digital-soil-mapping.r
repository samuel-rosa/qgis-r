options("repos"="http://cran.at.r-project.org/")
.libPaths("/home/alessandrorosa/Rlibs")
tryCatch(find.package("caret"), error = function(e) install.packages("caret", dependencies=TRUE))
library("caret")
tryCatch(find.package("snow"), error = function(e) install.packages("snow", dependencies=TRUE))
library("snow")
tryCatch(find.package("sf"), error = function(e) install.packages("sf", dependencies=TRUE))
library("sf")
tryCatch(find.package("raster"), error = function(e) install.packages("raster", dependencies=TRUE))
library("raster")
Observations <- st_read("/home/alessandrorosa/projects/software/qgis-r/data/vector/taxon-sample-val.shp", quiet = TRUE, stringsAsFactors = FALSE)
Response <- "taxon_sibc"
Weights <- "weights"
Validation <- "validation"
tempvar0 <- brick("/home/alessandrorosa/projects/software/qgis-r/data/raster/FT_29S54_.tif")
tempvar1 <- brick("/home/alessandrorosa/projects/software/qgis-r/data/raster/H3_29S54_.tif")
tempvar2 <- brick("/home/alessandrorosa/projects/software/qgis-r/data/raster/HN_29S54_.tif")
tempvar3 <- brick("/home/alessandrorosa/projects/software/qgis-r/data/raster/SN_29S54_.tif")
tempvar4 <- brick("/home/alessandrorosa/projects/software/qgis-r/data/raster/V3_29S54_.tif")
tempvar5 <- brick("/home/alessandrorosa/projects/software/qgis-r/data/raster/VN_29S54_.tif")
tempvar6 <- brick("/home/alessandrorosa/projects/software/qgis-r/data/raster/ZN_29S54_.tif")
tempvar7 <- brick("/home/alessandrorosa/projects/software/qgis-r/data/raster/geo50k.tif")
Covariates = c(tempvar0,tempvar1,tempvar2,tempvar3,tempvar4,tempvar5,tempvar6,tempvar7)
Model <- 4

# Load necessary packages ----
library(caret)
library(snow)

# Drop geometry and set response data type ----
Observations <- sf::st_drop_geometry(Observations)
if (is.character(Observations[[Response]])) {
  Observations[[Response]] <- as.factor(Observations[[Response]])
}

# Remove observations with NAs ----
na_idx <- complete.cases(Observations)
Observations <- Observations[na_idx, ]

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
} else if (any(Observations[[Validation]]) == 1) {
  validate <- TRUE
  idx <- which(Observations[[Validation]] == 1)
  val_data <- Observations[idx, ]
  Observations <- Observations[-idx, -which(colnames(Observations) == Validation)]
  n_val <- length(val_data[[Response]])
} else {
  validate <- FALSE
}

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

# rpart does not deal with parameters that are supposed to be passed to other models because they are checked
# against a list of valid function arguments.
# if (model == "rpart") {
model_fit <- caret::train(
  form = form,
  data = Observations,
  weights = Observations[[Weights]], 
  method = model,
  na.action = na.omit,
  tuneLength = 1,
  trControl = trainControl(method = "LOOCV")
)
# } else {
# model_fit <- caret::train(
# form = form, 
# data = Observations, 
# weights = Observations[[Weights]], 
# method = model, 
# na.action = na.omit,
# tuneLength = 1,
# trControl = trainControl(method = "LOOCV"),
# importance = TRUE, #rf
# prob.model = prob.model # svmRadial
# )
# }

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

writeRaster(Predictions, "/tmp/processing_e57569d77d1e4269a47daece2fd633b8/39f008ceef454cbca52c10c8b4056463/Predictions.tif", overwrite = TRUE)
writeRaster(Uncertainty, "/tmp/processing_e57569d77d1e4269a47daece2fd633b8/76fc91a871ea46de9f798ea57cb44c41/Uncertainty.tif", overwrite = TRUE)
write.csv(Metadata, "/tmp/processing_e57569d77d1e4269a47daece2fd633b8/2c58001c381846ae95810046d284dc16/Metadata.csv")
