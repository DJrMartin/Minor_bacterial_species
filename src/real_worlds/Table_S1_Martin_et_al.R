rm(list=ls())

### Packages.
packages <- readLines("requirements.txt")
invisible(lapply(packages, function(pkg){library(pkg, character.only = TRUE)}))

### data ================
load("data/data.rda")

## Normalisation ==================
X <- as.data.frame(t(apply(data$`Martin et al.`$count, 2, function(x) x/sum(x))))
# rowSums(X)

### Outcomes
FO <- data$`Martin et al.`$data$Lipides_SV1/9
Y <- as.factor(FO>0.4)

############################ PART 2 ############################
### data
x_clr <- as.matrix(compositions::clr(X))
x_PA <- X>0

source("functions/out_of_bag.R")
set.seed(2)
res.clr = out_of_bag_prediction(X, Y, transformed = TRUE, X_transformed = x_clr, cv = 10)
set.seed(2)
res.PA = out_of_bag_prediction(X, Y, transformed = TRUE, X_transformed = x_PA, cv = 10)
set.seed(2)
res.X = out_of_bag_prediction(as.matrix(X), Y, cv = 10)

df = rbind(paste0(round(apply(res.X$AUC.X, 2, function(x) mean(x[which(x!="NaN")])),2), " (±", round(apply(res.X$AUC.X, 2, function(x) sd(x[which(x!="NaN")])),2),")"),
      paste0(round(apply(res.X$AUC.R, 2, function(x) mean(x[which(x!="NaN")])),2), " (±", round(apply(res.X$AUC.R, 2, function(x) sd(x[which(x!="NaN")])),2),")"),
      paste0(round(apply(res.PA$AUC.X, 2, function(x) mean(x[which(x!="NaN")])),2), " (±", round(apply(res.PA$AUC.X, 2, function(x) sd(x[which(x!="NaN")])),2),")"),
      paste0(round(apply(res.PA$AUC.R, 2, function(x) mean(x[which(x!="NaN")])),2), " (±", round(apply(res.PA$AUC.R, 2, function(x) sd(x[which(x!="NaN")])),2),")"),
      paste0(round(apply(res.clr$AUC.X, 2, function(x) mean(x[which(x!="NaN")])),2), " (±", round(apply(res.clr$AUC.X, 2, function(x) sd(x[which(x!="NaN")])),2),")"),
      paste0(round(apply(res.clr$AUC.R, 2, function(x) mean(x[which(x!="NaN")])),2), " (±", round(apply(res.clr$AUC.R, 2, function(x) sd(x[which(x!="NaN")])),2),")"))

for(i in 1:dim(res.clr$AUC.R)[2]){
  print(t.test(res.X$AUC.X[,i], res.X$AUC.R[,i], paired = T)$p.value)
}

