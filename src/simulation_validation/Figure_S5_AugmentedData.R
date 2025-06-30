rm(list=ls())

## Library
library(yarrr)
library(caret)
library(randomForest)

## Data 
load("~/Desktop/for_github/data.rda")
load("~/Desktop/for_github/Simulation.RDATA")

NB_data.simulated = 40 # Nomber of simulated data for each method
X_ref <- data.frame(t(apply(data$PRJEB1220$count, 2, function(x) x/sum(x))))
Y_ref <- as.factor(data$PRJEB1220$data$disease)

n = nrow(X_ref)
p = ncol(X_ref)

AUC_from_BS = AUC_from_mida = AUC_from_ref = AUC_from_V1 = NULL
###. FIGURE S5: Augmented Data ###. 
set.seed(123)
for (iS in 1:NB_data.simulated){
  ## Reference model
  intraining <- createDataPartition(Y_ref, p=0.75, list=F)
  rf_train <- randomForest(Y_ref[intraining]~., X_ref[intraining, ])
  
  ## Version 1
  # colnames(Simulations$`Version 1`[[paste0("Xs_", iS)]]) = colnames(X_ref)
  # Y_sim <- predict(rf_train, Simulations$`Version 1`[[paste0("Xs_", iS)]])
  # rf_train_V1 <- randomForest(Y_sim~., Simulations$`Version 1`[[paste0("Xs_", iS)]])
  ## BiomeSampler
  Y_sim <- predict(rf_train, Simulations$BiomeSampler[[paste0("Xs_", iS)]])
  rf_train_BS <- randomForest(Y_sim~., Simulations$BiomeSampler[[paste0("Xs_", iS)]])
  ## MidaSim
  Y_sim <- predict(rf_train, Simulations$MiDAsim[[paste0("Xs_", iS)]])
  rf_train_mida <- randomForest(Y_sim~., Simulations$MiDAsim[[paste0("Xs_", iS)]])
  
  AUC_from_BS <- c(AUC_from_BS, pROC::auc(Y_ref[-intraining],predict(rf_train_BS, X_ref[-intraining, ], type = "prob")[,2]))
  AUC_from_mida <- c(AUC_from_mida, pROC::auc(Y_ref[-intraining],predict(rf_train_mida, X_ref[-intraining, ], type = "prob")[,2]))
  AUC_from_ref <- c(AUC_from_ref, pROC::auc(Y_ref[-intraining],predict(rf_train, X_ref[-intraining, ], type = "prob")[,2]))
  # AUC_from_V1 <- c(AUC_from_V1, pROC::auc(Y_ref[-intraining],predict(rf_train_V1, X_ref[-intraining, ], type = "prob")[,2]))
  print(iS)
  
}
boxplot(AUC_from_ref, AUC_from_BS, AUC_from_mida, col = c("grey",piratepal(palette = "pony")[2:3]), 
        axes=F, ylab = "AUC-ROC")
axis(1, at=1:3, labels = c("Ref.", "BiomeSampler", "MidaSim"))
axis(2)
segments(x0=2,x1=3, y0=0.95)
text(2.5, 0.98, "NS", xpd=NA)

