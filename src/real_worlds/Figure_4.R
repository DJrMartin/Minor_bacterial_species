rm(list=ls())
## IMPORTATIONS========================================
library(vegan)
library(ape)
library(dendextend)
library(SIBER)
library(cluster)
library(randomForest)

### data ================
load("data/data.rda")

## Normalisation ==================
X <- as.data.frame(t(apply(data$`Martin et al.`$count, 2, function(x) x/sum(x))))
# rowSums(X)

### Outcomes
FO <- data$`Martin et al.`$data$Lipides_SV1/9
Y <- as.factor(FO>0.4)

### Enterotype definition
bray_curtis <- vegan::vegdist(X) # Bray Curtis (dis)similarity
Z <- cutree(hclust(bray_curtis, method = "ward.D"), 2) # Enterotype

### Silhouette coef.
sil <- silhouette(Z, bray_curtis)
plot(sil)

### PCoA
res.pcoa = ape::pcoa(bray_curtis)
Coordinates <- res.pcoa$vectors
Eigenvalues <- res.pcoa$values[,3]*100

### Color definition
col.enterotype <- as.character(factor(Z, levels(as.factor(Z)), yarrr::piratepal(palette = "pony")[c(2,3)]))
Y.col = as.character(factor(Y, c("TRUE","FALSE"), c("#7F00FF90","#00640090")))

############################ PART 1 ############################
par(mar=c(4,4,2,2))
plot(Coordinates, xlab = "PCo1", ylab = "PCo2", col = Y.col, pch=16, axes = FALSE )

legend("bottomleft", legend=c("> 0.4 g/min" , "< 0.4 g/min"), fill=unique(Y_col), bty="n", cex=0.8)
title("PCoA - Rel. abundance data")

# Cluster drawing
for (i in unique(col.enterotype)){
  mu <- colMeans(Coordinates[which(col.enterotype==i),1:2])
  Sigma <- cov(Coordinates[which(col.enterotype==i),1:2]) 
  addEllipse(mu, Sigma, p.interval = 0.8, col = i, lty = 1, lwd=4)
}
axis(1)
axis(2)

# Random Forest pour avoir l'importance des species. 
rf <- randomForest(as.factor(Z) ~ . , data=X)

plot(log(colMeans(X)), log(apply(X,2, var)), 
     cex=rf$importance/2, axes=FALSE,
     xlim=c(-21, 5), ylim=c(-35, 5), main="Importance of each species on clustering",
     xlab = "Species' abundances (log scale)", ylab="Species' variability (log scale)")
axis(1)
axis(2)

############################ PART 2 ############################
### data
x_clr <- compositions::clr(X)
x_PA <- X>0

residuals.coda = residuals.clr = residuals.PA = clr = coda = PA = NULL
residuals.s.coda = residuals.s.clr = residuals.s.PA = clr.s = coda.s = PA.s = NULL
spe.coda = spe.PA = spe.r.PA = spe.r.coda = NULL

for(i in 1:50){
  intraining <- caret::createDataPartition(Y, p=0.8)$Resample1
  
  train = X[intraining,]
  test = X[-intraining,]
  bray_curtis <- vegan::vegdist(train) # Bray Curtis (dis)similarity
  Z <- as.factor(cutree(hclust(bray_curtis, method = "ward.D"), 2)) # Enterotype
  rf_ent <- randomForest(Z~., train)
  Z_train = Z
  Z_test = predict(rf_ent, test)
  
  R.coda = R.clr = R.PA = matrix(NA, nrow = nrow(x_clr), ncol = ncol(x_clr))
  for(j in levels(Z)){
    M_1 <- colMeans(X[intraining,][Z_train=="1",])
    M_2 <- colMeans(X[intraining,][Z_train=="2",])
    R.coda[intraining,][Z_train=="1",] = scale(X[intraining,][Z_train=="1",], center = M_1, scale = F)
    R.coda[intraining,][Z_train=="2",] = scale(X[intraining,][Z_train=="2",], center = M_2, scale = F)
    R.coda[-intraining,][Z_test=="1",] = scale(X[-intraining,][Z_test=="1",], center = M_1, scale = F)
    R.coda[-intraining,][Z_test=="2",] = scale(X[-intraining,][Z_test=="2",], center = M_2, scale = F)
  }
  
  #coda
  rf = randomForest::randomForest(Y[intraining]~., data = X[intraining,])
  coda <- c(coda, pROC::auc(Y[-intraining],predict(rf, X[-intraining,], type="prob")[,1], direction=">"))
  w <- which(rf$importance!=0)
  coda.s <- c(coda.s, summary(lm(log(rf$importance[w])~log(colMeans(X))[w]))$r.squared)
  spe.coda = cbind(spe.coda, rf$importance[,1])
  
  #clr
  # rf = randomForest::randomForest(Y[intraining]~., data=x_clr[intraining,])
  # clr <- c(clr, pROC::auc(Y[-intraining],predict(rf, x_clr[-intraining,], type="prob")[,1], direction=">"))
  # w <- which(rf$importance!=0)
  # clr.s <- c(clr.s, summary(lm(log(rf$importance[w])~log(colMeans(X))[w]))$r.squared)
  # 
  # #PA
  # rf = randomForest::randomForest(Y[intraining]~., data=x_PA[intraining,])
  # PA <- c(PA, pROC::auc(Y[-intraining],predict(rf, x_PA[-intraining,], type="prob")[,1], direction=">"))
  # w <- which(rf$importance!=0)
  # PA.s <- c(PA.s, summary(lm(log(rf$importance[w])~log(colMeans(X))[w]))$r.squared)
  # spe.PA = cbind(spe.PA, rf$importance[,1])
  
  #r
  rf = randomForest::randomForest(Y[intraining]~., data=R.coda[intraining,])
  residuals.coda <- c(residuals.coda, pROC::auc(Y[-intraining], predict(rf, R.coda[-intraining,], type="prob")[,1], direction=">"))
  w <- which(rf$importance!=0)
  residuals.s.coda <- c(residuals.s.coda, summary(lm(log(rf$importance[w])~log(colMeans(X))[w]))$r.squared)
  spe.r.coda = cbind(spe.r.coda, rf$importance[,1])
  
  # rf = randomForest::randomForest(Y[intraining]~., data=R.clr[intraining,])
  # residuals.clr <- c(residuals.clr, pROC::auc(Y[-intraining],predict(rf, R.clr[-intraining,], type="prob")[,1], direction=">"))
  # w <- which(rf$importance!=0)
  # residuals.s.clr <- c(residuals.s.clr, summary(lm(log(rf$importance[w])~log(colMeans(X))[w]))$r.squared)
  # 
  # rf = randomForest::randomForest(Y[intraining]~., data=R.PA[intraining,])
  # residuals.PA <- c(residuals.PA, pROC::auc(Y[-intraining],predict(rf, R.PA[-intraining,], type="prob")[,1], direction=">"))
  # w <- which(rf$importance!=0)
  # residuals.s.PA <- c(residuals.s.PA, summary(lm(log(rf$importance[w])~log(colMeans(X))[w]))$r.squared)
  # spe.r.PA = cbind(spe.r.PA, rf$importance[,1])
  
  
  print(i)
}

layout(matrix(c(1,2), nrow = 1))
par(mar=c(4,4,2,2))

color = rep(NA, 6)
color[seq(1,5,by=2)] = LaCroixColoR::lacroix_palette(type = "paired")[seq(2,6,by=2)]
color[seq(2,6,by=2)] = LaCroixColoR::lacroix_palette(type = "paired")[seq(1,5,by=2)]

boxplot(coda, residuals.coda, PA, residuals.PA, clr, residuals.clr, col=color, ylim=c(0,1.1), 
        axes=F, ylab="AUROC")
t.test(coda, residuals.coda)
abline(h=0.5, lty=2)
text(c(2,4,6),1.05, "a", cex=1)
axis(1, at=c(1.5,3.5,5.5),labels =rep('',3) )
text(c(1.5,3.5,5.5), y=rep(-0.22,4), c('Rel. Ab.',"Abs./pres.","CLR"), cex=0.8, srt=45, xpd=NA, pos=3)
axis(2, at=c(0.5,0.6, 0.7, 0.8, 0.9,1))


par(mar=c(4,5,2,2))
boxplot(coda.s,residuals.s.coda, PA.s,residuals.s.PA,clr.s,residuals.s.clr, col=color, ylim=c(0,0.7), 
        axes=F, ylab="Strenght of the relation between\nspecies importance and abundance")
axis(1, at=c(1.5,3.5,5.5),labels =rep('',3))
text(c(2,4,6),0.65, "a", cex=1)
text(c(1.5,3.5,5.5), y=rep(-0.13,4), c('Rel. Ab.',"Abs./pres.","CLR"), cex=0.8, srt=45, xpd=NA, pos=3)
axis(2)

plot.new()
legend("bottomleft", legend =rep("",6),fill=color, ncol=3, cex=0.8, bty="n", title="Type of data")
legend("bottomright", legend =c("Original data", "Guided transformed"), cex=0.8, bty = "n")
