##Modelling=group
##Spatial environmental modelling=name
##Observations=vector
##Target=field Observations
##Weights=field Observations
##Covariates=multiple raster
##Learner=string rf
##Predictions=output raster
##Uncertainty=output raster
##Metadata=output table

# Load necessary libraries ----
library(sp)
library(raster)
library(caret)
library(snow)

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

# Calibrate statistical learner ----
learner_fit <- train( 
  y = Observations[[Target]], weights = Observations[[Weights]],
  x = Observations@data[, !colnames(Observations@data) %in% c(Target, Weights)],
  method = Learner, tuneLength = 3,
  # preProcess = c("center", "scale"), 
  trControl= trainControl(method = "LOOCV"))

# Prepare for spatial predictions ----
if (is.numeric(Observations[[Target]])) {
  index <- 1
  type <- "raw"
} else {
  index <- 1:nlevels(Observations[[Target]])
  type <- "prob"
}

# Make spatial predictions ----
beginCluster()
prediction <- 
  clusterR(
    brick(Covariates), raster::predict, 
    args = list(model = learner_fit, type = type, index = index)
  )
endCluster()

# Compute predictions and prediction uncertainty ----
if (type == "prob") {
  Predictions <- as.factor(calc(x = prediction, fun = nnet::which.is.max))
  rat <- levels(Predictions)[[1]]
  rat$class <- levels(Observations[[Target]])
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
        paste("Predicted class (", paste(apply(rat, 1, paste, collapse = "="), collapse = "; "), ")", sep = "")),
      c("Uncertainty", 
        "Band 1 = Theoretical purity (0-1); Band 2 = Shannon entropy (0-1); Band 3 = Confusion index (0-1)"),
      c("Statistical learner", 
        paste(learner_fit$method[1], " = ", learner_fit$modelInfo$label[1], " (", learner_fit$modelType[1], ")", 
              sep = "")),
      c("Cross-validation", 
        paste("Accuracy = ", round(learner_fit$results$Accuracy[nrow(learner_fit$results)], 4), "; ",
              "Kappa = ", round(learner_fit$results$Kappa[nrow(learner_fit$results)], 4), 
              sep = ""))
    )
  colnames(Metadata) <- c("Item", "Description")
} else {
  Predictions <- prediction
  
  Uncertainty <- prediction # temporary
  
  Metadata <- 
    rbind(
      c("Predictions",  paste("Predicted values (", Target, ")", sep = "")),
      
      c("Uncertainty", paste("Predicted values (", Target, ")", sep = "")), # temporary
      
      c("Statistical learner",
      paste(learner_fit$method[1], " = ", learner_fit$modelInfo$label[1], " (", learner_fit$modelType[1], ")",
      sep = "")),
      c("Cross-validation",
      paste("RMSE = ", round(learner_fit$results$RMSE[nrow(learner_fit$results)], 4), "; ",
      "Rsquared = ", round(learner_fit$results$Rsquared[nrow(learner_fit$results)], 4), sep = ""))
    )
  colnames(Metadata) <- c("Item", "Description")
}

# Output ----
Predictions
Uncertainty
Metadata
