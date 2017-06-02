##Vector processing=group
##Convert data type to string=name
##Layer=vector
##Field=string
##Output=output vector

library(stringr)

Field <- stringr::str_split_fixed(string = Field, pattern = " ", n = Inf)
Field <- colnames(Layer@data) %in% Field
Layer@data[, Field] <- lapply(Layer@data[, Field], as.character)
Output <- Layer
Output
