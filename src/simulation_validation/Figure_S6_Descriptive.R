rm(list=ls())

## Library
library(yarrr)
## Data 
load("~/Desktop/for_github/data.rda")
load("~/Desktop/for_github/Simulation.RDATA")

NB_data.simulated = 40 # Nomber of simulated data for each method
X_ref <- data.frame(apply(data$PRJEB1220$count, 2, function(x) x/sum(x)))
p = nrow(X_ref)

###. FIGURE S4: Comparaison of the simulation methods. ###. 

# MEAN ==============================
layout(matrix(c(1:15), nrow=3, byrow = T))
par(mar=c(4,5,2,1))
RMSE=c()
for (m in 2:length(Simulations)){
  error = c()
  sim.covariates = c()
  for (iS in 1:NB_data.simulated){
    error = rbind(error, 
                  sqrt((rowMeans(X_ref)-colMeans(Simulations[[m]][[paste0("Xs_", iS)]]))^2))
    sim.covariates <- rbind(sim.covariates, colMeans(Simulations[[m]][[paste0("Xs_", iS)]]))
  }
  
  RMSE=cbind(RMSE, apply(error, 2, mean))
  plot(log(rowMeans(X_ref)), log(apply(sim.covariates, 2, mean)), 
       pch=16, col= piratepal(palette = "pony")[m], main=names(Simulations)[m], 
       xlab="Observed mean (log scale)", ylab="Simulated mean (log scale)\n(results from 40 simulations)")
  abline(0,1)
  print(paste("DONE!", m))
  
}
boxplot(log(RMSE), col= piratepal(palette = "pony")[2:5], axes=F, 
        ylab="RMSE (log scale)")
box()
axis(2)
# axis(1, at=1:5, labels=rep("", 5))
# text(x=1.5:5.5, rep(min(log(RMSE))*1.1, 5),pos=2,labels = names.method, srt=45, xpd=NA)

# PROP of ZEROs ==============================
# layout(matrix(c(1:6), nrow=1))
par(mar=c(4,5,2,1))
RMSE=c()
for (m in 2:length(Simulations)){
  error = c()
  sim.covariates = c()
  for (iS in 1:NB_data.simulated){
    obs <- rowSums(X_ref==0)/p
    simu <- colSums(Simulations[[m]][[paste0("Xs_", iS)]]==0)/p
    error = rbind(error, 
                  sqrt((obs-simu)^2))
    sim.covariates <- rbind(sim.covariates, simu)
  }
  
  RMSE=cbind(RMSE, apply(error, 2, mean))
  plot(log(obs), log(apply(sim.covariates, 2, mean)), 
       pch=16, col= piratepal(palette = "pony")[m],
       xlab="Observed 0 prop. (log scale)", ylab="Simulated 0 prop. (log scale)\n(results from 40 simulations)")
  abline(0,1)
  print(paste("DONE!", m))
  
}
wilcox.test(RMSE[,1], RMSE[,4])
boxplot(log(RMSE), col= piratepal(palette = "pony")[2:5], axes=F, 
        ylab="RMSE (log scale)")
box()
axis(2)
# axis(1, at=1:5, labels=rep("", 5))
# text(x=1.5:5.5, rep(min(log(RMSE))*1.1, 5),pos=2,labels = names.method, srt=45, xpd=NA)

# bray ==============================
# layout(matrix(c(1:6), nrow=1))
par(mar=c(4,5,2,1))
RMSE=c()
for (m in 2:length(Simulations)){
  error = c()
  sim.covariates = c()
  for (iS in 1:NB_data.simulated){
    obs <- as.numeric(vegan::vegdist(X_ref))
    simu <- as.numeric(vegan::vegdist(t(Simulations[[m]][[paste0("Xs_", iS)]])))
    error = rbind(error, 
                  sqrt((obs-simu)^2))
    sim.covariates <- rbind(sim.covariates, simu)
  }
  
  RMSE=cbind(RMSE, apply(error, 2, mean))
  plot(log(obs), log(apply(sim.covariates, 2, mean)), 
       pch=16, col= piratepal(palette = "pony")[m],
       xlab="Observed (dis)similarities (log scale)", ylab="Simulated (dis)similarities (log scale)\n(results from 40 simulations)")
  abline(0,1)
  print(paste("DONE!", m))
  
}
wilcox.test(RMSE[,2], RMSE[,3])
boxplot(log(RMSE), col= piratepal(palette = "pony")[2:5], axes=F, 
        ylab="RMSE (log scale)")
box()
axis(2)
