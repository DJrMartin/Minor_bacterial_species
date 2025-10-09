rm(list=ls())

### Packages.
packages <- readLines("requirements.txt")
invisible(lapply(packages, function(pkg){library(pkg, character.only = TRUE)}))
source("functions/out_of_bag.R")

## Data.
load("data/data.rda")

tss_data <- data.frame(t(apply(data$`Wu et al.`[,-c(1:4)], 1, function(x) x/sum(x)))) 
X <- tss_data[,colSums(tss_data)!=0]

id.study <- as.factor(data$`Wu et al.`$projectid)
Y <- as.factor(data$`Wu et al.`$group)

names.cohort.out = c("PRJEB33711", "PRJDB6133" , "PRJNA368966", "PRJNA431126", "PRJNA596546", "PRJNA398089")

res_R = res_X = NULL
for(out_of_bag in 1:length(names.cohort.out)){
  intraining <- which(id.study!=names.cohort.out[out_of_bag])

  y_train = Y[intraining]
  y_test = Y[-intraining]
  
  # Definition of the train and test sets for enterotypes
  train = X[intraining,]
  test = X[-intraining,]
  
  bray_curtis <- vegan::vegdist(train, "bray") 
  # Enterotype of the samples included in the train group
  Z_train <- as.factor(cutree(hclust(bray_curtis, method = "ward.D"), 3))
  
  # Prediction of the enterotypes from microbiome composition.
  rf_ent <- randomForest(Z_train ~ ., train)
  Z_test = predict(rf_ent, test)
  
  ## FUNCTION to compute the residuals of the out of bag
  # Initialize result matrices
  residuals.X <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
  
  # Loop over each group in enterotype
  for (grp in levels(Z_train)) {
    # Compute group-specific mean from training data
    M_grp <- colMeans(X[intraining,][Z_train == grp,])
    
    # Apply centering to training data
    residuals.X[intraining,][Z_train == grp,] <- scale(X[intraining,][Z_train == grp,], center = M_grp, scale = FALSE)
    # Apply centering to test data
    residuals.X[-intraining,][Z_test == grp,] <- scale(X[-intraining,][Z_test == grp,], center = M_grp, scale = FALSE)
  }
  
  set.seed(2)
  # PREDICTION
  # Prediction of Y from X
  rf_X = randomForest::randomForest(y_train~., data = X[intraining,], ntree = 1500, maxnodes = 5, mtry = 25)
  ## predictions
  pred.X.threshold <- predict(rf_X, X[-intraining,])
  pred.X <- predict(rf_X, X[-intraining,], type="prob")[,1]
  ## other metrics for X.
  metrics.X <- caret::confusionMatrix(pred.X.threshold, y_test,  mode = "everything")
  metrics.X <- c(metrics.X$overall[1], metrics.X$byClass[c(1,2,7)], cor(as.numeric(y_test), as.numeric(pred.X.threshold)))
  ## AUROC for X.
  res_X <- rbind(res_X, c(auc(y_test, pred.X), metrics.X))
  
  # Prediction on Y from Residuals
  rf_R = randomForest::randomForest(y_train~., data=residuals.X[intraining,], ntree = 1500, maxnodes = 5, mtry = 25)
  ## predictions
  pred.R <- predict(rf_R, residuals.X[-intraining,], type="prob")[,1]
  pred.R.threshold <- predict(rf_R, residuals.X[-intraining,])
  ## other metrics for r.
  metrics.R <- caret::confusionMatrix(pred.R.threshold, y_test, mode = "everything")
  metrics.R <- c(metrics.R$overall[1], metrics.R$byClass[c(1,2,7)], cor(as.numeric(y_test), as.numeric(pred.R.threshold)))
  ## AUROC for r.
  res_R <- rbind(res_R, round(c(auc(y_test, pred.R), metrics.R), 2))
}

res_X = round(res_X, 2)

for(i in 1:6){
  print(wilcox.test(res_X[,i], res_R[,i], paired = T)$p.value)
}

