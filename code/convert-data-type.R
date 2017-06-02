##Vector processing=group
##Convert data type=name
##Layer=vector
##Field=string
##Type=selectionReal;Integer;String
##Output=output vector

# Load necessary libraries ----
library(stringr)

# Identify fields ----
Field <- stringr::str_split_fixed(string = Field, pattern = " ", n = Inf)
if (!all(Field %in% colnames(Layer@data))) {
  error <- Field[which(!Field %in% colnames(Layer@data))]
  stop (paste("Attribute table does not containg Field", paste(error, collapse = " ")))
}
Field <- colnames(Layer@data) %in% Field

# Process fields ----
Type <- c("Real", "Integer", "String")[Type + 1]
if (Type == "Real") {
  Layer@data[, Field] <- sapply(Layer@data[, Field], as.numeric)
} else if (Type == "Integer") {
  Layer@data[, Field] <- sapply(Layer@data[, Field], as.integer)
} else if (Type == "String") {
  Layer@data[, Field] <- sapply(Layer@data[, Field], as.character)
}

# Output ----
Output <- Layer
Output
