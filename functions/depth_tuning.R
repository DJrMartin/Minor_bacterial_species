library(caret)
source("functions/out_of_bag.R")

depth_tuning <- function(X, Y, depth = 5, ntree = 1000, maxnodes = 20){
  # Initialisation
  # Use the model
  res <- out_of_bag_prediction(X, Y, cv = 1, ntree = ntree, maxnodes = maxnodes)
  # recorde results
  res.depth <- c(res$X, res$R)
  for(i in 2:depth){
    x = res$residual.data
    res <- out_of_bag_prediction(x, Y, cv = 1, ntree = ntree, maxnodes = maxnodes, D = "euclidean")
    # recorde results
    res.depth <- c(res.depth, res$R)
  }
  return(res.depth)
}