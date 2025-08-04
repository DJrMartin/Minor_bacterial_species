rm(list=ls())
######## TESTS
### data ================
load("data/data.rda")
source("functions/out_of_bag.R")

#==================================================================================
### MARTIN ET AL. (2025) ###
## X
X <- as.data.frame(t(apply(data$`Martin et al.`$count, 2, function(x) x/sum(x))))
# rowSums(X)
## Y
FO <- data$`Martin et al.`$data$Lipides_SV1/9
Y <- as.factor(FO>0.4)

res.martin = out_of_bag_prediction(X, Y, cv = 50)
# boxplot(res.martin$X, res.martin$R) ; t.test(res.martin$X, res.martin$R)
#==================================================================================