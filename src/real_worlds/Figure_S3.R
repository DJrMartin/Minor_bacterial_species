rm(list=ls())
library(SIBER)
library(ape)
library(tidyverse)
library(vegan)
library(dendextend)
library(cluster)
library(mclust)

### data ================
load("data/data.rda")

### Normalisation ================
X <- data$PRJEB39223$count
Y <- data$PRJEB39223$data$BMI

### Preprocess the tree.
labels=tidyr::separate(data.frame(S=data$`Martin et al.`$tree$tip.label), 'S', c("CAG",'K', 'P', 'C', "O", "F", "G", "S"), sep = "__" )
data$`Martin et al.`$tree$tip.label=labels$S
rownames(labels)=labels$S
filt_tree <- ape::keep.tip(data$`Martin et al.`$tree, intersect(data$`Martin et al.`$tree$tip.label, rownames(X)))

### Color of the outcomes variable.
colfunc <- colorRampPalette(c('purple',"ivory3",'darkgreen'))(6)
Y_col = rep(NA, length(Y))
for (i in 1:6){
  Y_col[which(Y>c(10,18,25,30,35,50)[i] & Y<=c(10,18,25,30,35,50)[i+1])]=colfunc[i]
}

### DATA FRAME ================
x_bray <- (vegan::vegdist(t(X)))
x_clr <- (dist(compositions::clr(t(X))))
x_unifrac <- (picante::unifrac(t(X), filt_tree))
x_PA <- (vegan::vegdist(t(X>0)))

### Guided transformation on CLR-transformed data
GRP_clr <- hclust(vegan::vegdist(t(X)), method = "ward.D")
Z1 <- cutree(GRP_clr, 2)
mod = lm(compositions::clr(t(X))~as.factor(Z1))
R = data.frame(mod$residuals)
x_R <- (dist(R))

DF <- list(x_bray, x_clr, x_PA, x_unifrac ,x_R)
names.method <- c("Bray-Curtis", "CLR", "Pres/abs.", "Unifrac", "Guided Transformation")

##########################################################
################ FIGURE S3 ###############################
##########################################################

layout(matrix(c(1,1,1,2,3,3,3,4,5,5,5,6,
                1,1,1,2,3,3,3,4,5,5,5,6,
                1,1,1,2,3,3,3,4,5,5,5,6,
                7,7,7,8,9,9,9,10,11,11,11,12,
                7,7,7,8,9,9,9,10,11,11,11,12,
                7,7,7,8,9,9,9,10,11,11,11,12), ncol=12, byrow = T))
for(i in 1:5){
  par(mar=c(4,6,3,1))
  x.y <- ape::pcoa(DF[[i]])$vectors
  Z2 <- cutree(hclust(DF[[i]], method = "ward.D"), 2)
  plot(x.y, col=as.character(Y_col), pch=16, ylim=c(min(x.y[,2]), max(x.y[,2])*1.2),
       xlab="PCo1", ylab="PCo2", main=names.method[i])
  ## cluster 1
  mu <- colMeans(x.y[which(Z2==1),1:2])
  Sigma <- cov(x.y[which(Z2==1),1:2]) 
  addEllipse(mu, Sigma, p.interval = 0.5, col = "red", lty = 1, lwd=3)
  ## cluster 2
  mu <- colMeans(x.y[which(Z2==2),1:2])
  Sigma <- cov(x.y[which(Z2==2),1:2]) 
  addEllipse(mu, Sigma, p.interval = 0.5, col = "blue", lty = 1, lwd=3)
  ## Boxplot
  par(mar=c(4,4,3,0))
  boxplot(Y[which(Z2==1)], Y[which(Z2==2)], col=c("red", 'blue'), ylim=c(18, 50), ylab="BMI", axes=F)
  axis(2)
  text(1.5, 48,paste(round(t.test(Y[which(Z2==1)], Y[which(Z2==2)])$p.value, 3)), cex=1)
  segments(x0= 1, x1=2, y0=46)
}

plot.new()
legend("top", title =expression("Body Mass Index"),legend=c("]10-18]", "]18-25]","]25-30]","]30-35]","]35-50]"), ncol=3, fill=colfunc, bty="n")
legend("bottom", title =expression(paste("Estimated cluster ",Z)), legend=c("1", '2'), fill=c("red", "blue"), bty="n")

