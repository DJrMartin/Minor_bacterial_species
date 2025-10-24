rm(list=ls())

### Packages.
packages <- readLines("requirements.txt")
invisible(lapply(packages, function(pkg){library(pkg, character.only = TRUE)}))

### FOR DISTINCT PREPROCESS
# Preprocess =========================
load("data/data.rda")

x = data.frame(t(apply(data$PRJEB1220$count, 2, function(x) x/sum(x))))
# rowSums(X)

x.clr = data.frame(compositions::clr(x)) # CLR transformation
x.rclr = data.frame(decostand(x, method = "rclr"))
colnames(x.rclr) = colnames(x.clr)
x.pa = data.frame(x>0)

Y <- factor(data$PRJEB1220$data$disease, c("Health","Irritable Bowel Syndrome"),c("HLT", "IBD"))

x$Class = Y
x.clr$Class = Y
x.rclr$Class = Y
x.pa$Class = Y

NAMES <- c("Rel. Ab.", "clr", "rclr", "Pres.abs")
DF = list(x, x.clr, x.rclr, x.pa)

# Conditions =========================
cv_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = multiClassSummary)

empty.list = c()

res = list(list(empty.list, empty.list), 
           list(empty.list, empty.list),
           list(empty.list, empty.list),
           list(empty.list, empty.list))

for(i in 1:length(DF)){
  X = DF[[i]]
  AUC = res.interpretability = c()
  
  set.seed(12)
  for(j in 1:10){
    
    # train
    train_index <- createDataPartition(Y, p = 0.4, list = FALSE)
    train_data <- X[train_index, ]
    # test
    test_index <- sample(as.matrix(1:396)[-train_index,], 100)
    test_data <- X[test_index,]
    
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
    AUC = rbind(AUC, c(pROC::auc(test_data$Class ,predict(models$GBM, test_data, "prob")[,2]), 
                                 pROC::auc(test_data$Class ,predict(models$RF, test_data, "prob")[,2]),
                                 pROC::auc(test_data$Class ,predict(models$LASSO, test_data, "prob")[,2]),
                                 pROC::auc(test_data$Class ,predict(models$PLS, test_data, "prob")[,2])))
    
    # Variable Importance
    importance_list <- lapply(models, varImp)
    interpretability = c()
    
    for(m in 1:4){
      if(i == 4){
        I <- gsub("TRUE", "", x = names(sort(as.matrix(importance_list[[m]]$importance)[,1], decreasing=T)[1:20])) 
      }else{I <- names(sort(as.matrix(importance_list[[m]]$importance)[,1], decreasing=T)[1:20])}
      interpretability = cbind(interpretability, colMeans(x[,I]))
    }
    
    res.interpretability = rbind(res.interpretability, interpretability)
    
  }
  res[[i]][[1]] = AUC
  res[[i]][[2]] = res.interpretability
  
  # Avancment of the analyses.
  print(paste("done", i))
}

color = LaCroixColoR::lacroix_palette(type = "paired")[c(2, 4, 6, 8, 12)]

layout(matrix(c(1,2,3), nrow=1))
par(mar=c(5,5,4,2))
at.bx <- c(0.6, 0.8, 1, 1.2)

### FIGURE I
boxplot(res[[1]][[1]], axes=F, col=color[1], boxwex=0.4, xlim = c(0.5, 4.5),
        ylim=c(0.5,1), ylab="AUC-ROC", main="Predictive performance", at = at.bx[1]:3.6)
for(i in 2:4){
  boxplot(res[[i]][[1]], axes = F, col = color[i],
          cex.main = 0.8, boxwex = 0.4,add=T, at = at.bx[i]:4.2)
}
axis(1, 1:4,labels=F)
axis(2)
text(1:4, y=0.43, c("GBM","RF", "LASSO", "PLS"), xpd=NA, srt=60)
legend("bottomleft", legend=c("Rel. ab.","CLR", "r-CLR","Abs./pres."), fill = color, 
       bty='n', cex=1)

### FIGURE J
boxplot(res[[1]][[2]], axes=F, col=color[1], boxwex=0.4, xlim = c(0.5, 4.5),
        ylim=c(0.5,1), ylab="Abundances (log scale) of the\n20 most important species", 
        main="Important species defined by models", at = at.bx[1]:3.6)
for(i in 2:4){
  boxplot(res[[i]][[2]], axes = F, col = color[i],
          cex.main = 0.8, boxwex = 0.4,add=T, at = at.bx[i]:4.2)
}
axis(1, 1:4, labels=F)
axis(2)
text(1:4, y=-17, c("GBM","RF", "LASSO", "PLS"), xpd=NA, srt=60)
legend("bottomleft", legend=c("Rel. ab.","CLR", "r-CLR","Abs./pres."), fill = color, 
       bty='n', cex=1)

### FIGURE K
rf = randomForest::randomForest(Class ~ ., data = x)
plot(log(colMeans(as.matrix(for.next))), log(rf$importance), cex=0.5, xlab="Relative Abundance of species\n(log-scale)", ylab="Importance of the species\nin the RF model",
     col=color[1], pch = 3,ylim=c(-8,2), main="Relation between species abundance and\ntheir importance in the RF model.")
w = which(rf$importance!=0)
summary(lm(log(colMeans(as.matrix(for.next)))[w]~log(rf$importance)[w]))
abline(coef(lm(log(rf$importance)[w]~log(colMeans(as.matrix(for.next)))[w])))
legend("bottom", legend=c(expression(paste("r = 0.31 ; ", r^2," = 0.56, p < 0.001"))),cex=1, bty="n")
