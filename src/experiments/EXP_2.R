rm(list=ls())
# DATA =========================
library(compositions)
library(randomForest)
library(SIBER)
library(yarrr)
library(FactoMineR)

load("data/Simulation.RDATA")
load("data/data.rda")
# Preprocess =========================
X = data.frame(t(data$PRJEB1220$count))
BC <- vegan::vegdist(X) # BC
X.clr <- clr(X) # CLR transformation

# # Enterotype estimation =========================
# # Z_1 <- as.numeric(cutree(hclust(BC, method = "ward.D"), 2)) 
Z_for_R <- as.numeric(cutree(hclust(BC, method = "ward.D"), 3))
# # Residuals estimation =========================
R <- lm(as.matrix(X.clr)~as.factor(Z_for_R))$residuals # Calcul of the residuals.
# Host physiology =========================
Y <- as.factor(read.csv('~/microbiota/host_physiology.csv', row.names = 1)[,1])
# Selection of the Y's support. =========================
J <- which(apply(X, 2, var) > quantile(apply(X, 2, var), p=0.4))
# Random Forest to compute Y' =========================
# rf_1 = randomForest(as.factor(Z_1)~., X[,J])
rf_2 = randomForest(as.factor(Y)~., R[,-J])

########################
####### PARTIE 2 ####### 
########################
N = 299
set.seed(1)
res.coda = res.residuals = res.clr = res.PA  = NULL
res.inter.coda = res.inter.clr = res.inter.r = res.inter.PA = NULL
for(nb in 1:49){
  # nb=2
  # (A) Preprocess of the X' =========================
  X.sim <- Simulations$BiomeSampler[[nb]]
  X.sim[which(X.sim<0)]=0
  X.sim = as.data.frame(t(apply(X.sim, 1, function(x) x/sum(x)*100)))
  colnames(X.sim) <- colnames(X)
  
  # (B) Enterotype of X' =========================
  Z_1 <- cutree(hclust(vegan::vegdist(X.sim), method = "ward.D2"), 2) # Clustering
  prerequis <- length(which(Z_1==1))<length(which(Z_1==2))
  if(prerequis==TRUE){
    Z_1.corrected = Z_1
    Z_1[which(Z_1.corrected==1)]=2 ; Z_1[which(Z_1.corrected==2)]=1
  }
  Z_1_R <- cutree(hclust(vegan::vegdist(X.sim), method = "ward.D2"), 3) # Clustering
  # plot(ape::pcoa(vegan::vegdist(X.sim))$vectors, col=as.factor(Z_1_R))
  
  # (C) Take a subsequent of the data.
  ech.train <- caret::createDataPartition(Z_1_R, p=N/400, list=F)
  # TRAIN ==============
  entero.train = Z_1[ech.train]
  coda.train <- X.sim[ech.train,] # Compositionnal data.
  PA.train <- X.sim[ech.train,]>0
  clr.train <- data.frame(clr(coda.train)) # CLR
  Z_1_R <- cutree(hclust(vegan::vegdist(coda.train), method = "ward.D2"), 3) # Clustering
  r.pred.train <- data.frame(lm(as.matrix(clr.train)~as.factor(Z_1_R))$residuals)
  r.train <- data.frame(lm(as.matrix(PA.train)~as.factor(Z_1_R))$residuals)
  
  # TEST ==============
  ech.test <- as.numeric(names(sample(Z_1[-ech.train], 100)))+1
  entero.test = Z_1[ech.test]
  coda.test <- X.sim[ech.test,] # Compositionnal data.
  PA.test <- X.sim[ech.test,]>0
  clr.test <- data.frame(clr(coda.test)) # CLR
  Z1_test <- cutree(hclust(vegan::vegdist(coda.test), method = "ward.D2"), 3) # Clustering
  r.pred.test <- data.frame(lm(as.matrix(clr.test)~as.factor(Z1_test))$residuals)
  r.test <- data.frame(lm(as.matrix(PA.test)~as.factor(Z1_test))$residuals) # Residuals
  
  # plot(ape::pcoa(vegan::vegdist(coda.test))$vectors, col=as.factor(Z1_test))
  # plot(ape::pcoa(vegan::vegdist(coda.train))$vectors, col=as.factor(Z_1_R))
  
  ### REST OF THE ALGORITHMS
  p.coda = p.r = p.clr = p.PA  = NULL
  spe.coda = spe.r = spe.clr = spe.PA = NULL
  for(w_1 in c(0, 0.08, 0.15)){
    # w_1 =0
    # (D) Simulation of Y_{sim}
    # a.Probability to belong to the Prevotella enterotype
    # p_1 <- as.numeric(predict(rf_1, rbind(coda.train, coda.test), type = "prob")[,2])
    p_1 <- c(as.numeric(entero.train-1), as.numeric(entero.test)-1)
    # b.Probability to belong to the IBD class.
    p_2 <- as.numeric(predict(rf_2, rbind(r.pred.train, r.pred.test), type = "prob")[,2])
    # c. Y Weighted Average
    w_2 = 1-w_1
    Y_sim <- (w_1*(p_1) +  w_2*(p_2)) > 0.5

    Y_sim.train = as.factor(Y_sim[1:dim(r.train)[1]])
    Y_sim.test = as.factor(Y_sim[(dim(r.train)[1]+1):(dim(r.train)[1]+100)])
    
    # (E) randomforest
    rf_r <- randomForest(Y_sim.train~., r.train)
    rf_bray <- randomForest(Y_sim.train~., coda.train)
    rf_clr <- randomForest(Y_sim.train~., clr.train)
    rf_PA <- randomForest(Y_sim.train~., PA.train)

    # (G) interpretability
    # (G) interpretability
    spe.coda = c(spe.coda, mean(as.numeric(log(apply(coda.train, 2, var))[names(sort(rf_bray$importance[,1], decreasing = T))[1:20]])))
    spe.clr = c(spe.clr, mean(as.numeric(log(colMeans(coda.train))[names(sort(rf_clr$importance[,1], decreasing = T))[1:20]])))
    spe.r = c(spe.r, mean(as.numeric(log(colMeans(coda.train))[names(sort(rf_r$importance[,1], decreasing = T))[1:20]])))
    spe.PA = c(spe.PA, mean(as.numeric(log(colMeans(coda.train))[names(sort(rf_PA$importance[,1], decreasing = T))[1:20]])))
    
    # (F) Results of predictive performance on test
    p.r = c(p.r, pROC::auc(Y_sim.test, predict(rf_r, r.test, type='prob')[,2]))
    p.coda = c(p.coda, pROC::auc(Y_sim.test, predict(rf_bray, coda.test, type='prob')[,2]))
    p.clr = c(p.clr, pROC::auc(Y_sim.test, predict(rf_clr, clr.test, type='prob')[,2]))
    p.PA = c(p.PA, pROC::auc(Y_sim.test, predict(rf_PA, PA.test, type='prob')[,2]))
    
  }
  res.coda <- rbind(res.coda, c(p.coda))
  res.residuals <- rbind(res.residuals, c(p.r))
  res.clr <- rbind(res.clr, c(p.clr))
  res.PA <- rbind(res.PA, c(p.PA))
  
  res.inter.coda = rbind(res.inter.coda, c(spe.coda))
  res.inter.clr = rbind(res.inter.clr, c(spe.clr))
  res.inter.r = rbind(res.inter.r, c(spe.r))
  res.inter.PA = rbind(res.inter.PA, c(spe.PA))

  print(paste(nb, "done ! "))
}

layout(matrix(c(1,2), nrow = 1))
par(mar=c(4,5,2,2))
color = LaCroixColoR::lacroix_palette(type = "paired")[seq(2,8, by=2)]
color[4] = LaCroixColoR::lacroix_palette(type = "paired")[3]

##### FIGURE 1
DF=list(res.coda, res.PA, res.clr , res.residuals)
plot(1,1,xlim=c(0.5, 3.5),ylim=c(0.5, 1), pch=15, col='white', 
     axes=F, xlab="", ylab='AUC-ROC')
axis(2)
axis(1, at=1:3, labels=rep("", 3))
text(1:3, y = 0.4 ,c("Low", "Medium", "High"), xpd=NA)

positions = c(0.8, 0.9, 1, 1.1, 1.2)
for (i in 1:4){
  M = apply(DF[[i]], 2, mean)
  SD = apply(DF[[i]], 2, function(x) sd(x))
  points(x = positions[i]:(positions+2)[i] ,M, ylim=c(-50,0), pch=15, col = color[i], type = "b")
  segments(x0 = positions[i]:(positions+2)[i], y0 = M-SD , y1 = M+SD, col = color[i])
}


# Noms des groupes pour les comparaisons
group_names <- c("CODA", "PA","CLR", "R")

# Pour chaque colonne (1 à 3)
for (col in 1:3) {
  cat("-----\nColonne", col, "\n")

  # Extraire les données pour cette colonne
  values <- unlist(lapply(DF, function(x) x[, col]))
  groups <- factor(rep(group_names, each = nrow(DF[[1]])))

  # Test de Kruskal-Wallis
  kw_result <- kruskal.test(values ~ groups)
  cat("Kruskal-Wallis p-value:", kw_result$p.value, "\n")

  if (kw_result$p.value < 0.05) {
    cat("=> p < 0.05, test de Wilcoxon par paire :\n")

    # Comparaisons par paire
    pairwise_result <- pairwise.wilcox.test(values, groups, p.adjust.method = "BH")
    print(pairwise_result)
  } else {
    cat("=> p >= 0.05, pas de test pairwise effectué.\n")
  }
}
# 300 et 150
text(c(1), y = 0.98, c("d"), cex = 0.8)
# text(c(1:3), y = 0.98, c("b,c',d'", "c,d'", "d"), cex = 0.8)
# text(c(1:3), y = 0.98, c("c,d'", "c,d'", "d"), cex = 0.8)

##### FIGURE 2
DF=list(res.inter.coda, res.inter.PA,res.inter.clr, res.inter.r)
plot(1,1,xlim=c(0.5, 3.5),ylim=c(-6,0), pch=15, col='white', 
     axes=F, xlab="", ylab='Mean (log scale) of the signal\ndriving the clustering')
axis(2)
axis(1, at=1:3, labels=rep("", 3))
text(1:3, y = -7 ,c("Low", "Medium", "High"), xpd=NA)

positions = c(0.8, 0.9, 1, 1.1, 1.2)
for (i in 1:4){
  M = apply(DF[[i]], 2, mean)
  SD = apply(DF[[i]], 2, function(x) sd(x))
  points(x = positions[i]:(positions+2)[i] ,M, pch=15, col=color[i], type = "b")
  segments(x0 = positions[i]:(positions+2)[i], y0 = M-SD , y1 = M+SD, col=color[i])
}

