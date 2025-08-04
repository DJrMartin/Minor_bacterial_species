### data ================
load("data/data.rda")
source("functions/depth_tuning.R")

#===============================================================================
### WU ET AL. (2024) ###
X <- as.data.frame(t(apply(data$`Wu et al.`[,-c(1:4)], 1, function(x) x/sum(x))))
# rowSums(X)
## Y
Y <- as.factor(data$`Wu et al.`$group)

res.wu = depth_tuning(X, Y)
#===============================================================================