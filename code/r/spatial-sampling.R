##Point patterns=group
##Spatial sampling=name
##Layer=vector
##Type=selection Random;Regular;Stratified;Nonaligned;Hexagonal;Clustered
##Size=number
##Clusters=number
##Output=output vector

# Identify sample type
Type <- Type + 1
type <- c("random", "regular", "stratified", "nonaligned", "hexagonal", "clustered")
type <- type[Type]

# Sample
Output <- sp::spsample(x = Layer, n = Size, type = type, iter = 10, nclusters = Clusters)
Output <- sp::SpatialPointsDataFrame(Output, as.data.frame(Output))
