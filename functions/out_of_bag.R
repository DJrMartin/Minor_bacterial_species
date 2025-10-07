library(randomForest)
library(pROC)

out_of_bag_prediction <- function(X, Y, transformed = FALSE, X_transformed = NULL, cv = 20, ntree = 1500, 
                                  maxnodes = 5, mtry = 5, distance = "bray", k = 2, p = 0.7){
  res_R = res_X = NULL
  res.imp.R = res.imp.X = NULL
  # cross validation
  for(cv in 1:cv){
    intraining <- caret::createDataPartition(Y, p = p, list = F)
    # Target variable
    y_train = Y[intraining]
    y_test = Y[-intraining]
    
    # Definition of the train and test sets for enterotypes
    train.ent = X[intraining,]
    test.ent = X[-intraining,]
    
    bray_curtis <- vegan::vegdist(train.ent, distance) 
    # Enterotype of the samples included in the train group
    Z_train <- as.factor(cutree(hclust(bray_curtis, method = "ward.D"), k))
    
    # Prediction of the enterotypes from microbiome composition.
    rf_ent <- randomForest(Z_train ~ ., train.ent)
    Z_test = predict(rf_ent, test.ent)
    
    if (transformed==TRUE){
      X_transformed = X_transformed
    }else{
      X_transformed = X
    }
    
    ## FUNCTION to compute the residuals of the out of bag
    # Initialize result matrices
    residuals.X <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
    
    # Loop over each group in enterotype
    for (grp in levels(Z_train)) {
      # Compute group-specific mean from training data
      M_grp <- colMeans(X_transformed[intraining,][Z_train == grp,])
      
      # Apply centering to training data
      residuals.X[intraining,][Z_train == grp,] <- scale(X_transformed[intraining,][Z_train == grp,], center = M_grp, scale = FALSE)
      # Apply centering to test data
      residuals.X[-intraining,][Z_test == grp,] <- scale(X_transformed[-intraining,][Z_test == grp,], center = M_grp, scale = FALSE)
    }
    
    # PREDICTION
    # Prediction of Y from X
    rf_X = randomForest::randomForest(y_train~., data = X_transformed[intraining,], ntree = ntree, maxnodes = maxnodes, mtry = mtry)
    res_X <- c(res_X, auc(y_test, predict(rf_X, X_transformed[-intraining,], type="prob")[,1], direction=">"))
    # Prediction on Y from Residuals
    rf_R = randomForest::randomForest(y_train~., data=residuals.X[intraining,], ntree = ntree, maxnodes = maxnodes, mtry = mtry)
    res_R <- c(res_R, auc(y_test, predict(rf_R, residuals.X[-intraining,], type="prob")[,1], direction=">"))
  
    # IMPORTANCE
    # Importance in X or X_transformed
    w <- which(rf_X$importance!=0)
    res.imp.X <- c(res.imp.X, summary(lm(log(rf_X$importance[w])~log(colMeans(X))[w]))$r.squared)
    # Importance in residuals
    w <- which(rf_R$importance!=0)
    res.imp.R <- c(res.imp.R, summary(lm(log(rf_R$importance[w])~log(colMeans(X))[w]))$r.squared)
  
  }
  
  # return the resutls.
  return(list(AUC.X = res_X, AUC.R = res_R, Imp.X = res.imp.X, Imp.R = res.imp.R))
}
