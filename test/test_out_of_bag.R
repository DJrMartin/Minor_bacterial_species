rm(list=ls())
######## TESTS
### data ================
load("data/data.rda")
source("functions/out_of_bag.R")

#==================================================================================
### MARTIN ET AL. (2025) ###
## X
X <- as.data.frame(t(apply(data$`Martin et al.`$count, 2, function(x) x/sum(x))))
x_clr <- as.matrix(compositions::clr(X))
x_PA <- X>0
# rowSums(X)
## Y
FO <- data$`Martin et al.`$data$Lipides_SV1/9
Y <- as.factor(FO>0.4)

set.seed(2)
res = out_of_bag_prediction(X, Y, transformed = FALSE, cv = 5)
#==================================================================================
