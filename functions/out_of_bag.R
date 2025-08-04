library(randomForest)
library(pROC)

out_of_bag_prediction <- function(X, Y, cv = 30, ntree = 1500, maxnodes = 5, D = "bray", k = 2){
  res_R = res_X = NULL
  # cross validation
  for(cv in 1:cv){
    intraining <- caret::createDataPartition(Y, p = 0.8, list = F)
    # Definition of the train and test sets
    train = X[intraining,]
    test = X[-intraining,]
   
    y_train = Y[intraining]
    y_test = Y[-intraining]
    if(D == "bray"){# Bray Curtis (dis)similarity
      bray_curtis <- vegan::vegdist(train) 
      # Enterotype of the samples included in the train group
      Z_train <- as.factor(cutree(hclust(bray_curtis, method = "ward.D"), k))
    }
    if(D != "bray"){
      bray_curtis <- dist(train) 
      # Enterotype of the samples included in the train group
      Z_train <- as.factor(cutree(hclust(bray_curtis, method = "ward.D"), k))
    }
    
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
    
    # Prediction of Y from X
    rf_X = randomForest::randomForest(y_train~., data = train, ntree = ntree, maxnodes = maxnodes)
    res_X <- c(res_X, auc(y_test, predict(rf_X, test, type="prob")[,1], direction=">"))
    
    # Prediction on Y from Residuals
    rf_R = randomForest::randomForest(y_train~., data=residuals.X[intraining,])
    res_R <- c(res_R, auc(y_test, predict(rf_R, residuals.X[-intraining,], type="prob")[,1], direction=">"))
  
  }
  
  # return the resutls.
  return(list(X = res_X, R = res_R, residual.data = residuals.X))
}
