rm(list=ls())

### Packages.
packages <- readLines("requirements")
invisible(lapply(packages, function(pkg){library(pkg, character.only = TRUE)}))

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

legend("bottomleft", legend=c("> 0.4 g/min" , "< 0.4 g/min"), fill = unique(Y.col), bty="n", cex=0.8)
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
x_clr <- as.matrix(compositions::clr(X))
x_PA <- X>0

source("functions/out_of_bag.R")
set.seed(2)
res.clr = out_of_bag_prediction(X, Y, transformed = TRUE, X_transformed = x_clr, cv = 20)
set.seed(2)
res.PA = out_of_bag_prediction(X, Y, transformed = TRUE, X_transformed = x_PA, cv = 20)
set.seed(2)
res.X = out_of_bag_prediction(as.matrix(X), Y, cv = 20)

layout(matrix(c(1,2), nrow = 1))
par(mar=c(4,4,2,2))

color = rep(NA, 6)
color[seq(1,5,by=2)] = LaCroixColoR::lacroix_palette(type = "paired")[seq(2,6,by=2)]
color[seq(2,6,by=2)] = LaCroixColoR::lacroix_palette(type = "paired")[seq(1,5,by=2)]

boxplot(res.X$AUC.X, res.X$AUC.R, res.PA$AUC.X, res.PA$AUC.R, res.clr$AUC.X, res.clr$AUC.R, col=color, ylim=c(0,1.1), 
        axes=F, ylab="AUROC")
t.test(res.X$AUC.X, res.X$AUC.R, paired = T) ; t.test(res.PA$AUC.X, res.PA$AUC.R, paired = T) ; t.test(res.clr$AUC.X, res.clr$AUC.R, paired = T)

abline(h=0.5, lty=2)
text(c(2,4,6),1.05, "a", cex=1)
axis(1, at=c(1.5,3.5,5.5),labels =rep('',3) )
text(c(1.5,3.5,5.5), y=rep(-0.22,4), c('Rel. Ab.',"Abs./pres.","CLR"), cex=0.8, srt=45, xpd=NA, pos=3)
axis(2, at=c(0.5,0.6, 0.7, 0.8, 0.9,1))

par(mar=c(4,5,2,2))
boxplot(res.X$Imp.X, res.X$Imp.R, res.PA$Imp.X, res.PA$Imp.R, res.clr$Imp.X, res.clr$Imp.R, col=color, ylim=c(0,0.7), 
        axes=F, ylab="Strenght of the relation between\nspecies importance and abundance")
t.test(res.X$Imp.X, res.X$Imp.R, paired = T) ; t.test(res.PA$Imp.X, res.PA$Imp.R, paired = T) ; t.test(res.clr$Imp.X, res.clr$Imp.R, paired = T)

axis(1, at=c(1.5,3.5,5.5),labels =rep('',3))
text(c(2,4,6),0.65, "a", cex=1)
text(c(1.5,3.5,5.5), y=rep(-0.13,4), c('Rel. Ab.',"Abs./pres.","CLR"), cex=0.8, srt=45, xpd=NA, pos=3)
axis(2)

plot.new()
legend("bottomleft", legend =rep("",6),fill=color, ncol=3, cex=0.8, bty="n", title="Type of data")
legend("bottomright", legend =c("Original data", "Guided transformed"), cex=0.8, bty = "n")
