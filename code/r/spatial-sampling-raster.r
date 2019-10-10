##Point patterns=group
##Spatial sampling (raster)=name
##Layer=raster
##Type=selection Random;Regular;Stratified;Nonaligned;Hexagonal;Clustered
##Size=number
##Clusters=number
##Seed=number
##Output=output vector

# Identify sample type
Type <- Type + 1
type <- c("random", "regular", "stratified", "nonaligned", "hexagonal", "clustered")
type <- type[Type]

# Convert raster to SpatialPixelsDataFrame
Layer <- as(Layer, 'SpatialPixelsDataFrame')

# Sample
set.seed(Seed)
Output <- sp::spsample(x = as(Layer, 'SpatialPixelsDataFrame'), n = Size, type = type, iter = 10, nclusters = Clusters)
Output <- sp::SpatialPointsDataFrame(Output, as.data.frame(Output))
