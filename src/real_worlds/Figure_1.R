rm(list=ls())
library(caret)
library(gbm)
library(randomForest)
library(e1071)
library("MLmetrics")
library(compositions)
library(yarrr)

### FOR DISTINCT PREPROCESS
# Preprocess =========================
load("data/data.rda")
X = data.frame(t(apply(data$PRJEB1220$count, 2, function(x) x/sum(x))))
# rowSums(X)
X.clr <- data.frame(clr(X)) # CLR transformation
Y <- factor(data$PRJEB1220$data$disease, c("Health","Irritable Bowel Syndrome"),c("HLT", "IBD"))

for.next = X
R <- data.frame(apply(X>0, 2, as.numeric))
# R <- data.frame(lm(as.matrix(X)~as.factor(cutree(hclust(vegan::vegdist(X), method="ward.D"),2)))$residuals)
X$Class = Y
X.clr$Class = Y
R$Class = Y

# Conditions =========================
cv_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = multiClassSummary)

AUC.coda = res.coda = c()
AUC.clr = res.clr = c()
AUC.R = res.R = c()

set.seed(12)
for(i in 1:20){
  train_index <- createDataPartition(Y, p = 0.4, list = FALSE)
  train_data <- X[train_index, ]
  test_index <- sample(as.matrix(1:396)[-train_index,], 100)
  # intersect(train_index, test_index)
  test_data <- X[test_index,]
  
  ###### CODA ====================================================================
  # Model training with hyperparameter tuning for Coda.
  models <- list()
  # GBM
  models$GBM <- train(Class ~ ., data = train_data, method = "gbm",
                      trControl = cv_control, verbose=FALSE)
  # Random Forest
  models$RF <- train(Class ~ ., data = train_data, method = "rf", 
                     trControl = cv_control)
  # Lasso Regression (glmnet with alpha = 1)
  models$LASSO <- train(Class ~ ., data = train_data, method = "glmnet", 
                        trControl = cv_control, tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 1, length = 10)))
  # pls Regression 
  models$PLS <- train(Class ~ ., data = train_data, method = "pls", trControl = cv_control)
  # Model Comparison
  AUC.coda = rbind(AUC.coda, c(pROC::auc(test_data$Class ,predict(models$GBM, test_data, "prob")[,2]), 
                               pROC::auc(test_data$Class ,predict(models$RF, test_data, "prob")[,2]),
                               pROC::auc(test_data$Class ,predict(models$LASSO, test_data, "prob")[,2]),
                               pROC::auc(test_data$Class ,predict(models$PLS, test_data, "prob")[,2])))
  
  # Variable Importance
  importance_list <- lapply(models, varImp)
  interpretability=c()
  
  for(m in 1:4){
    I <- names(sort(as.matrix(importance_list[[m]]$importance)[,1], decreasing=T)[1:20])
    interpretability = cbind(interpretability , log(colMeans(X[,I])))
  }
  
  res.coda = rbind(res.coda, interpretability)
  
  ###### R
  train_data <- R[train_index, ]
  test_data <- R[test_index,]
  # Model training with hyperparameter tuning for Coda.
  models <- list()
  # GBM
  models$GBM <- train(Class ~ ., data = train_data, method = "gbm",
                      trControl = cv_control, verbose=FALSE)
  # Random Forest
  models$RF <- train(Class ~ ., data = train_data, method = "rf", 
                     trControl = cv_control)
  # Lasso Regression (glmnet with alpha = 1)
  models$LASSO <- train(Class ~ ., data = train_data, method = "glmnet", 
                        trControl = cv_control, tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 1, length = 10)))
  # pls Regression 
  models$PLS <- train(Class ~ ., data = train_data, method = "pls", trControl = cv_control)
  # Model Comparison
  AUC.R = rbind(AUC.R, c(pROC::auc(test_data$Class ,predict(models$GBM, test_data, "prob")[,2]), 
                         pROC::auc(test_data$Class ,predict(models$RF, test_data, "prob")[,2]),
                         pROC::auc(test_data$Class ,predict(models$LASSO, test_data, "prob")[,2]),
                         pROC::auc(test_data$Class ,predict(models$PLS, test_data, "prob")[,2])))
  
  # Variable Importance
  importance_list <- lapply(models, varImp)
  interpretability=c()
  
  for(m in 1:4){
    I <- names(sort(as.matrix(importance_list[[m]]$importance)[,1], decreasing=T)[1:20])
    interpretability = cbind(interpretability , log(colMeans(X[,I])))
  }
  
  res.R = rbind(res.R, interpretability)
  
  # ====================================================================================
  train_data <- X.clr[train_index, ]
  test_data <- X.clr[test_index, ]
  # Model training with hyperparameter tuning for CLR
  models <- list()
  # GBM
  models$GBM <- train(Class ~ ., data = train_data, method = "gbm",
                      trControl = cv_control, verbose=FALSE)
  # Random Forest
  models$RF <- train(Class ~ ., data = train_data, method = "rf", 
                     trControl = cv_control)
  # Lasso Regression (glmnet with alpha = 1)
  models$LASSO <- train(Class ~ ., data = train_data, method = "glmnet", 
                        trControl = cv_control, tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 1, length = 10)))
  # pls Regression 
  models$PLS <- train(Class ~ ., data = train_data, method = "pls", trControl = cv_control)
  
  AUC.clr = rbind(AUC.clr, c(pROC::auc(test_data$Class ,predict(models$GBM, test_data, "prob")[,2]), 
                             pROC::auc(test_data$Class ,predict(models$RF, test_data, "prob")[,2]),
                             pROC::auc(test_data$Class ,predict(models$LASSO, test_data, "prob")[,2]),
                             pROC::auc(test_data$Class ,predict(models$PLS, test_data, "prob")[,2])))
  
  # Variable Importance
  importance_list <- lapply(models, varImp)
  interpretability=c()
  
  for(m in 1:4){
    I <- names(sort(as.matrix(importance_list[[m]]$importance)[,1], decreasing=T)[1:20])
    interpretability = cbind(interpretability , log(colMeans(X[,I])))
  }
  
  res.clr = rbind(res.clr, interpretability)
  print(paste("done", i))
}

layout(matrix(c(1,2,3), nrow=1))
par(mar=c(5,5,4,2))
color = LaCroixColoR::lacroix_palette(type = "paired")[seq(2,8, by=2)]
color[5] = LaCroixColoR::lacroix_palette(type = "paired")[7]

boxplot(AUC.coda, axes=F, col=color[1], boxwex=0.4,
        ylim=c(0.5,1), ylab="AUC-ROC", main="Predictive performance",at = 0.6:3.6)
boxplot(AUC.R, axes=F, col=color[2],
        cex.main=0.8, boxwex=0.4,add=T,at = 0.8:3.8)
boxplot(AUC.clr, axes=F, col=color[3],
        cex.main=0.8, boxwex=0.4,add=T,at = 1:4)

axis(1, 1:4,labels=F)
axis(2)
text(1:4, y=0.43, c("GBM","RF", "LASSO", "PLS"), xpd=NA, srt=60)
legend("bottomleft", legend=c("Rel. ab.","Abs./pres.","CLR"), fill = color, 
       bty='n', cex=1)

boxplot(res.coda, axes=F, col=color[1], boxwex=0.4,
        ylab="Abundances (log scale) of the\n20 most important species", main="Important species defined by models",at = 0.6:3.6)
boxplot(res.R, axes=F, col=color[2],
        cex.main=0.8, boxwex=0.4,add=T,at = 0.8:3.8)
boxplot(res.clr, axes=F, col=color[3],
        cex.main=0.8, boxwex=0.4,add=T,at = 1:4)

axis(1, 1:4, labels=F)
axis(2)
text(1:4, y=-17, c("GBM","RF", "LASSO", "PLS"), xpd=NA, srt=60)
legend("bottomleft", legend=c("Rel. Ab.","Abs./pres.","CLR"), fill = color, 
       bty='n', cex=1)

rf = randomForest::randomForest(Class ~ ., data = X)

plot(log(colMeans(as.matrix(for.next))), log(rf$importance), cex=0.5, xlab="Relative Abundance of species\n(log-scale)", ylab="Importance of the species\nin the RF model",
     col=color[1], pch = 3,ylim=c(-8,2), main="Relation between species abundance and\ntheir importance in the RF model.")
w = which(rf$importance!=0)
summary(lm(log(colMeans(as.matrix(for.next)))[w]~log(rf$importance)[w]))
abline(coef(lm(log(rf$importance)[w]~log(colMeans(as.matrix(for.next)))[w])))
legend("bottom", legend=c(expression(paste("r = 0.31 ; ", r^2," = 0.56, p < 0.001"))),cex=1, bty="n")
