rm(list=ls())

### Packages.
packages <- readLines("requirements.txt")
invisible(lapply(packages, function(pkg){library(pkg, character.only = TRUE)}))

load("data/Simulation.RDATA")
load("data/data.rda")

# Preprocess =========================
X = data.frame(t(data$PRJEB1220$count))
BC <- vegan::vegdist(X) # BC
X.clr <- clr(X) # CLR transformation

# Enterotype estimation =========================
Z_for_R <- as.numeric(cutree(hclust(BC, method = "ward.D"), 3))
# # Residuals estimation =========================
R <- lm(as.matrix(X.clr)~as.factor(Z_for_R))$residuals # Calcul of the residuals.
# Host physiology =========================
Y <- as.factor(data$PRJEB1220$data$disease)
# Selection of the Y's support. =========================
J <- which(apply(X, 2, var) > quantile(apply(X, 2, var), p=0.4))
# Random Forest to compute Y' =========================
rf_2 = randomForest(as.factor(Y)~., R[,-J])

########################
####### PARTIE 2 ####### 
########################
n = c(75, 150, 299)
for(N in n){
  set.seed(1)
  res.coda = res.residuals = res.clr = res.rclr = res.PA  = NULL
  res.inter.coda = res.inter.clr = res.inter.rclr = res.inter.r = res.inter.PA = NULL
  for(nb in 1:40){
    # nb=2
    # (A) Preprocess of the X' =========================
    X.sim <- Simulations$BiomeSampler[[nb]]
    X.sim[which(X.sim<0)]=0
    X.sim = as.data.frame(t(apply(X.sim, 1, function(x) x/sum(x)*100)))
    colnames(X.sim) <- colnames(X)
    
    # (B) Enterotype of X' =========================
    Z_1 <- cutree(hclust(vegan::vegdist(X.sim), method = "ward.D"), 2) # Clustering
    prerequis <- length(which(Z_1==1))<length(which(Z_1==2))
    if(prerequis==TRUE){
      Z_1.corrected = Z_1
      Z_1[which(Z_1.corrected==1)]=2 ; Z_1[which(Z_1.corrected==2)]=1
    }
    Z_1_R <- cutree(hclust(vegan::vegdist(X.sim), method = "ward.D"), 3) # Clustering
    # plot(ape::pcoa(vegan::vegdist(X.sim))$vectors, col=as.factor(Z_1_R))
    
    # (C) Take a subsequent of the data.
    ech.train <- caret::createDataPartition(Z_1_R, p=N/400, list=F)
    
    # TRAIN ==============
    entero.train = Z_1[ech.train]
    
    coda.train <- X.sim[ech.train,] # Compositionnal data.
    PA.train <- X.sim[ech.train,]>0
    clr.train <- data.frame(clr(coda.train)) # CLR
    rclr.train <- data.frame(decostand(coda.train, "rclr")) # r-CLR
    
    Z1_train <- cutree(hclust(vegan::vegdist(coda.train), method = "ward.D"), 3) # Clustering
    
    r.pred.train <- data.frame(lm(as.matrix(clr.train)~as.factor(Z1_train))$residuals)
    r.train <- data.frame(lm(as.matrix(PA.train)~as.factor(Z1_train))$residuals)
    
    # TEST ==============
    ech.test <- as.numeric(names(sample(Z_1[-ech.train], 100)))+1
    entero.test = Z_1[ech.test]
    
    coda.test <- X.sim[ech.test,] # Compositionnal data.
    PA.test <- X.sim[ech.test,]>0
    clr.test <- data.frame(clr(coda.test)) # CLR
    rclr.test <- data.frame(decostand(coda.test, "rclr")) # r-CLR
    
    Z1_test <- cutree(hclust(vegan::vegdist(coda.test), method = "ward.D"), 3) # Clustering
    
    r.pred.test <- data.frame(lm(as.matrix(clr.test)~as.factor(Z1_test))$residuals)
    r.test <- data.frame(lm(as.matrix(PA.test)~as.factor(Z1_test))$residuals) # Residuals
    
    ### REST OF THE ALGORITHMS
    p.coda = p.r = p.clr = p.PA = p.rclr = NULL
    spe.coda = spe.r = spe.clr = spe.PA = spe.rclr = NULL
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
      rf_rclr <- randomForest(Y_sim.train~., rclr.train)
      rf_PA <- randomForest(Y_sim.train~., PA.train)
      
      # (G) interpretability
      spe.coda = c(spe.coda, mean(as.numeric(log(apply(coda.train, 2, var))[names(sort(rf_bray$importance[,1], decreasing = T))[1:20]])))
      spe.clr = c(spe.clr, mean(as.numeric(log(colMeans(coda.train))[names(sort(rf_clr$importance[,1], decreasing = T))[1:20]])))
      spe.rclr = c(spe.rclr, mean(as.numeric(log(colMeans(coda.train))[names(sort(rf_rclr$importance[,1], decreasing = T))[1:20]])))
      spe.r = c(spe.r, mean(as.numeric(log(colMeans(coda.train))[names(sort(rf_r$importance[,1], decreasing = T))[1:20]])))
      spe.PA = c(spe.PA, mean(as.numeric(log(colMeans(coda.train))[names(sort(rf_PA$importance[,1], decreasing = T))[1:20]])))
      
      # (F) Results of predictive performance on test
      p.r = c(p.r, pROC::auc(Y_sim.test, predict(rf_r, r.test, type='prob')[,2]))
      p.coda = c(p.coda, pROC::auc(Y_sim.test, predict(rf_bray, coda.test, type='prob')[,2]))
      p.clr = c(p.clr, pROC::auc(Y_sim.test, predict(rf_clr, clr.test, type='prob')[,2]))
      p.rclr = c(p.rclr, pROC::auc(Y_sim.test, predict(rf_rclr, rclr.test, type='prob')[,2]))
      p.PA = c(p.PA, pROC::auc(Y_sim.test, predict(rf_PA, PA.test, type='prob')[,2]))
      
    }
    
    res.coda <- rbind(res.coda, c(p.coda))
    res.residuals <- rbind(res.residuals, c(p.r))
    res.clr <- rbind(res.clr, c(p.clr))
    res.rclr <- rbind(res.rclr, c(p.rclr))
    res.PA <- rbind(res.PA, c(p.PA))
    
    res.inter.coda = rbind(res.inter.coda, c(spe.coda))
    res.inter.clr = rbind(res.inter.clr, c(spe.clr))
    res.inter.rclr = rbind(res.inter.rclr, c(spe.rclr))
    res.inter.r = rbind(res.inter.r, c(spe.r))
    res.inter.PA = rbind(res.inter.PA, c(spe.PA))
    
    print(paste(nb, "done ! "))
  }
  
  DF_1 = list(res.coda, res.clr, res.rclr, res.PA, res.residuals)
  DF_2 = list(res.inter.coda, res.inter.clr, res.inter.rclr, res.inter.PA, res.inter.r)
  res = list(DF_1, DF_2)
  
  # save(res, file = paste0("data/figure_3_exp2_", N, ".rda"))
}

# Statistique ============================================ 
# Noms des groupes pour les comparaisons
group_names <- c("CODA","CLR", "r_clr", "PA", "R")

# load(file = paste0("data/figure_3_exp2_299.rda"))
# Pour chaque colonne (1 à 3)
for (col in 1:3) {
  cat("-----\nColonne", col, "\n")
  
  # Extraire les données pour cette colonne
  values <- unlist(lapply(res[[1]], function(x) x[, col]))
  groups <- factor(rep(group_names, each = nrow(res[[1]][[1]])))
  
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

n = c(75, 150, 299)
for(N in n){
  # load(file = paste0("data/figure_3_exp2_", N, ".rda"))
  pdf(file = paste0("figures/figure_3_exp2_", N, ".pdf"), width = 6, height = 3)
 
  layout(matrix(c(1,2), nrow=1))
  par(mar=c(3,5,2,2))
  
  color = LaCroixColoR::lacroix_palette(type = "paired")[c(2, 4, 6, 12)]
  color[5] = LaCroixColoR::lacroix_palette(type = "paired")[11]
  
  ##### FIGURE 1
  plot(1,1,xlim=c(0.5, 3.5),ylim=c(0.5, 1), pch=15, col='white', 
       axes=F, xlab="", ylab='AUC-ROC')
  axis(2)
  axis(1, at=1:3, labels=rep("", 3))
  text(1:3, y = 0.4 ,c("Low", "Medium", "High"), xpd=NA)
  
  positions = c(0.8, 0.9, 1, 1.1, 1.2)
  for (i in 1:5){
    M = apply(res[[1]][[i]], 2, mean)
    SD = apply(res[[1]][[i]], 2, function(x) sd(x))
    points(x = positions[i]:(positions+2)[i] ,M, ylim=c(-50,0), pch=15, col = color[i], type = "b")
    segments(x0 = positions[i]:(positions+2)[i], y0 = M-SD , y1 = M+SD, col = color[i])
  }
  if(N==75){text(c(1:3), y = 0.98, c("c,d'", "c,d'", "d'"), cex = 0.8)}
  if(N==150){text(c(1:2), y = 0.98, c("c,d'", "c,d'"), cex = 0.8)}
  if(N==299){text(c(1), y = 0.98, c("d"), cex = 0.8)}
   
  ##### FIGURE 2
  plot(1,1, xlim = c(0.5, 3.5), ylim = c(-6,0), pch = 15, col = 'white', 
       axes = F, xlab = "", ylab = 'Mean (log scale) of the signal\ndriving the clustering')
  axis(2)
  axis(1, at = 1:3, labels = rep("", 3))
  text(1:3, y = -7 , c("Low", "Medium", "High"), xpd = NA)
  
  positions = c(0.8, 0.9, 1, 1.1, 1.2)
  for (i in 1:5){
    M = apply(res[[2]][[i]], 2, mean)
    SD = apply(res[[2]][[i]], 2, function(x) sd(x))
    points(x = positions[i]:(positions+2)[i] ,M, pch=15, col=color[i], type = "b")
    segments(x0 = positions[i]:(positions+2)[i], y0 = M-SD , y1 = M+SD, col=color[i])
  }
  dev.off()
}


