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

### Preprocess the tree.
labels=tidyr::separate(data.frame(S=data$`Martin et al.`$tree$tip.label), 'S', c("CAG",'K', 'P', 'C', "O", "F", "G", "S"), sep = "__" )
data$`Martin et al.`$tree$tip.label=labels$S
rownames(labels)=labels$S
filt_tree <- ape::keep.tip(data$`Martin et al.`$tree, intersect(data$`Martin et al.`$tree$tip.label, colnames(X)))

########################
####### PARTIE 1 ####### 
########################
# Enterotype estimation =========================
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
N = 75
res.coda = res.residuals = res.clr = res.PA = res.unifrac = NULL
res.inter.coda = res.inter.clr = res.inter.uni = res.inter.r = res.inter.PA = NULL

for(nb in 1:50){
  # nb=1
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
  
  # (C) Take a subsequent of the data.
  ech.train <- caret::createDataPartition(Z_1, p=N/400, list=F)
  # TRAIN ==============
  entero.train = Z_1[ech.train]
  # table(entero.train)
  coda.train <- X.sim[ech.train,] # Compositionnal data.
  P.A.train <- X.sim[ech.train,]>0 # Compositionnal data.
  # dim(train)
  clr.train <- data.frame(clr(coda.train)) # CLR
  # rclr.train <- data.frame(t(apply(coda.train, 1, function(x) transform(x, coda.train)))) # robust CLR
  Z_1_R <- cutree(hclust(vegan::vegdist(coda.train), method = "ward.D2"), 3) # Clustering
  r.train <- data.frame(lm(as.matrix(P.A.train)~as.factor(Z_1_R))$residuals) # Residuals
  r.pred.train <- data.frame(lm(as.matrix(clr.train)~as.factor(Z_1_R))$residuals) # Residuals
  d.Unifrac <- picante::unifrac(coda.train, filt_tree) ## Distance UniFrac.

  p.coda = p.r = p.clr = p.PA = p.unifrac = NULL
  spe.coda = spe.r = spe.clr = spe.PA = spe.unifrac = NULL
  
  for(w_1 in c(0, 0.08, 0.15)){
    # (D) Simulation of Y_{sim}
    # a.Probability to belong to the Prevotella enterotype
    # p_1 <- as.numeric(predict(rf_1, coda.train, type = "prob")[,2])
    p_1 <- as.numeric(entero.train-1)
    # b.Probability to belong to the IBD class.
    p_2 <- as.numeric(predict(rf_2, r.pred.train, type = "prob")[,2])
    # c. Y Weighted Average
    w_2 = 1-w_1
    Y_sim <- (w_1*(p_1) +  w_2*(p_2)) > 0.5
    # plot(ape::pcoa(vegan::vegdist(coda.train))$vectors, col=as.factor(entero.train), pch=as.numeric(Y_sim))

    # (E) Clustering of coda, clr, rclr and r.
    k_r <- kmeans(dist(r.train),2)$cluster
    k_clr <- kmeans(dist(clr.train),2)$cluster
    k_PA <- kmeans(dist(P.A.train),2)$cluster
    k_unifrac <- kmeans(dist(d.Unifrac),2)$cluster
    
    # (F) if Y_{sim} = P_2.
    # p.value between entero and Y.
    p.coda = c(p.coda, log(chisq.test(table(data.frame(Y = Y_sim, P_2 = p_1>0.5)))$p.value))
    # p.value between clustering on r and Y weigthed.
    p.r = c(p.r, log(chisq.test(table(data.frame(Y = Y_sim, k = k_r)))$p.value))
    # p.value between clustering on clr and Y weigthed.
    p.clr = c(p.clr, log(chisq.test(table(data.frame(Y = Y_sim, k = k_clr)))$p.value))
    # p.value between clustering on rclr and Y weigthed.
    p.PA = c(p.PA, log(chisq.test(table(data.frame(Y = Y_sim, k = k_PA)))$p.value))
    # p.value between clustering on Unifrac and Y weigthed.
    p.unifrac = c(p.unifrac, log(chisq.test(table(data.frame(Y = Y_sim, k = k_unifrac)))$p.value))

    # (G) interpretability
    rf <- randomForest::randomForest(as.factor(entero.train)~.,data.frame(coda.train))
    spe.coda = c(spe.coda, mean(as.numeric(log(colMeans(coda.train))[names(sort(rf$importance[,1], decreasing = T))[1:20]])))
    
    rf <- randomForest::randomForest(as.factor(k_clr)~.,data.frame(coda.train))
    spe.clr = c(spe.clr, mean(as.numeric(log(colMeans(coda.train))[names(sort(rf$importance[,1], decreasing = T))[1:20]])))
    
    rf <- randomForest::randomForest(as.factor(k_unifrac)~.,data.frame(coda.train))
    spe.unifrac = c(spe.unifrac, mean(as.numeric(log(colMeans(coda.train))[names(sort(rf$importance[,1], decreasing = T))[1:20]])))
    
    rf <- randomForest::randomForest(as.factor(k_PA)~.,data.frame(coda.train))
    spe.PA = c(spe.PA, mean(as.numeric(log(colMeans(coda.train))[names(sort(rf$importance[,1], decreasing = T))[1:20]])))
    
    rf <- randomForest::randomForest(as.factor(k_r)~.,data.frame(coda.train))
    spe.r = c(spe.r, mean(as.numeric(log(colMeans(coda.train))[names(sort(rf$importance[,1], decreasing = T))[1:20]])))
  }
  res.coda <- rbind(res.coda, c(p.coda))
  res.residuals <- rbind(res.residuals, c(p.r))
  res.clr <- rbind(res.clr, c(p.clr))
  res.PA <- rbind(res.PA, c(p.PA))
  res.unifrac <- rbind(res.unifrac, c(p.unifrac))
  
  res.inter.coda = rbind(res.inter.coda, c(spe.coda))
  res.inter.clr = rbind(res.inter.clr, c(spe.clr))
  res.inter.uni = rbind(res.inter.uni, c(spe.unifrac))
  res.inter.r = rbind(res.inter.r, c(spe.r))
  res.inter.PA = rbind(res.inter.PA, c(spe.PA))

  print(paste(nb, "done ! "))
}

layout(matrix(c(1,2), nrow=1))
color = LaCroixColoR::lacroix_palette(type = "paired")[seq(2,8, by=2)]
color[5] = LaCroixColoR::lacroix_palette(type = "paired")[3]
par(mar=c(2,5,2,2))
##### FIGURE 1
DF=list(res.coda, res.PA, res.clr, res.unifrac , res.residuals)
plot(1,1,xlim=c(0.5, 3.5),ylim=c(-30,0), pch=15, col='white', 
     axes=F, xlab="", ylab='p value (log scale)')
axis(2)
axis(1, at=1:3, labels=rep("", 3))
text(1:3, y = -35 ,c("Low", "Medium", "High"), xpd=NA)

positions = c(0.8, 0.9, 1, 1.1, 1.2)
for (i in 1:5){
  M = apply(DF[[i]], 2, mean)
  SD = apply(DF[[i]], 2, function(x) sd(x))
  points(x = positions[i]:(positions+2)[i] ,M, ylim=c(-50,0), pch=15, col=color[i], type = "b")
  segments(x0 = positions[i]:(positions+2)[i], y0 = M-SD , y1 = M+SD, col=color[i])
}
abline(h=log(0.05), lty=2)

##### FIGURE 2
DF=list(res.inter.coda, res.inter.PA, res.inter.clr, res.inter.uni, res.inter.r)
plot(1,1,xlim=c(0.5, 3.5),ylim=c(-4,0), pch=15, col='white', 
     axes=F, xlab="", ylab='Mean (log scale) of the signal\ndriving the clustering')
axis(2)
axis(1, at=1:3, labels=rep("", 3))
text(1:3, y = -4.6 ,c("Low", "Medium", "High"), xpd=NA)

positions = c(0.8, 0.9, 1, 1.1, 1.2)
for (i in 1:5){
  M = apply(DF[[i]], 2, mean)
  SD = apply(DF[[i]], 2, function(x) sd(x))
  points(x = positions[i]:(positions+2)[i] ,M, pch=15, col=color[i], type = "b")
  segments(x0 = positions[i]:(positions+2)[i], y0 = M-SD , y1 = M+SD, col=color[i])
}
plot.new()
legend("bottomleft", legend=c("Rel. Ab.", "Abs./Pres.", "CLR", "Unifrac", "Guided Transf."), fill=color[1:5],
       bty="n", cex=0.8, ncol=5)
