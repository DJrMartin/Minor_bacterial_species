rm(list=ls())

### Packages.
packages <- readLines("requirements.txt")
invisible(lapply(packages, function(pkg){library(pkg, character.only = TRUE)}))
source(file = "functions/composition.R")

load("data/data.rda")

X <- data.frame(t(apply(data$PRJEB1220$count, 2, function(x) x/sum(x))))
rowSums(X)

## Distance de Bray Curtis
BC <- vegan::vegdist(X)
coord = ape::pcoa(BC)$vectors[,1:10]

k=3
latent_k.BC <- as.factor(cutree(hclust(BC, method = "ward.D"), k))
latent_k.col <- as.character(factor(latent_k.BC, levels(as.factor(latent_k.BC)), piratepal(palette = "pony")[c(1:k)]))

Y <- as.factor(data$PRJEB1220$data$disease)
Y_col = as.character(factor(Y, c("Irritable Bowel Syndrome", "Health"), c("#7F00FF70","#00640070")))

pdf(file = "figures/figure_1_x.pdf", width = 12.5, height = 3.5)
layout(matrix(c(1,1,2,3,3,4,4,
                1,1,2,3,3,4,4), nrow=2, byrow = T))
### First part of the figure, describing the enterotypes and the PCoA. ============================================
par(mar=c(4,4,2,2))
# Projection de l'ensemble des études.
plot(coord, xlab="PCo1", ylab="PCo2", col=Y_col, pch=16, axes = FALSE )
legend("bottomright", legend=c("IBD" , "HLT"), fill=unique(Y_col), bty="n", cex=0.8)

title("PCoA - Rel. abundance data")
for (i in unique(latent_k.col)){
  mu <- colMeans(coord[which(latent_k.col==i),1:2])
  Sigma <- cov(coord[which(latent_k.col==i),1:2]) 
  addEllipse(mu, Sigma, p.interval = 0.8, col = i, lty = 1, lwd=4)
}
axis(1)
axis(2)

### Second part of the figure, describing the composition. ============================================
#### Computed the composition.
res.genus <- genus.compo(X, data$`Martin et al.`$tree)[[1]]

res.per.group = c()
for(i in 1:k){
  res.per.group <- cbind(res.per.group, colMeans(res.genus[which(latent_k.BC==i),]))
}

w.b <- which(apply(res.per.group, 1, var)>0.0001) # select the most expressed genus.
length(w.b)

composition.genus = res.per.group[w.b,]
Others = 1 - colSums(res.per.group[w.b,])
composition.genus = rbind(composition.genus, Others)

#### The barplot.
par(mar=c(4,4,2,2))
# dt %>% plot(horiz = T, main='', axes=F)
plot(x=rep(1, 3), y = seq(1,k,by=1),  col=yarrr::piratepal("pony")[1:3], 
     axes=F, xlab="", cex=5, pch=15, ylim=c(0.5,3.5), ylab="Enterotypes")
for(i in 1:k){text(x=1, y=i, xpd=NA,c("Mix-ET","P-ET","B-ET")[i])}
for(i in 1:3){text(x=1, y=i-0.35, xpd=NA,
                         paste(table(Y[which(latent_k.BC==i)])[1],"/", table(Y[which(latent_k.BC==i)])[2]))}
text(1, 3.5, "Healthy / IBD", xpd=NA)
text(1, 3.35, "p-value = 0.91", xpd=NA)
chisq.test(Y, latent_k.col)

par(mar=c(4,4,2,2))
barplot(composition.genus, horiz = T, col = yarrr::piratepal("info2"), 
        main='Gut microbiota composition',  border=NA)

### Thrid part of the figure, describing the composition. ============================================
rf <- randomForest(latent_k.BC~., data = X)
plot(log(colMeans(X)), log(apply(X, 2, var)), 
     cex=rf$importance/2, axes=FALSE,
     xlim=c(-21, 5), ylim=c(-35, 5), main="Importance of each species on clustering",
     xlab = "Species' abundances (log scale)", ylab="Species' variability (log scale)")
axis(1)
axis(2)
dev.off()

# Fourth part of the figure, describing the most important species. ============================================ 
# varImpPlot(rf, n=10)

# ALTERNATIVE CLUSTERING =================================================================================
k = 2 # number of alternative clustering
## transformed data
x.clr = dist(compositions::clr(X))
x.unifrac = picante::unifrac(X, genus.compo(X, data$`Martin et al.`$tree)[[2]])
x.rclr = vegan::vegdist(X, method = "robust.aitchison")
x.pa = vegan::vegdist(X>0, method = "jaccard")

NAMES <- c("clr", "Unweighted Unifrac", "rclr", "Pres.abs")
DF = list(x.clr, x.unifrac, x.rclr, x.pa)

for(i in 1:length(DF)){
  pdf(file = paste0("figures/figure_1_", NAMES[i], ".pdf"), width = 12.5, height = 3.5)
  layout(matrix(c(1,1,2,3,3,4,4,
                  1,1,2,3,3,4,4), nrow=2, byrow = T))
  
  # PART 1 ================================================
  par(mar=c(4,4,2,2))
  X.Y = ape::pcoa(DF[[i]])$vectors
  latent_k <- as.factor(cutree(hclust(DF[[i]], method = "ward.D2"), 2))
  
  Z_for_correction = latent_k
  if(diff(table(data.frame(latent_k, Y))[,2])<0){latent_k[Z_for_correction==1] = 2 ; latent_k[Z_for_correction==2] = 1}
  latent_k.col <- as.character(factor(latent_k, levels(as.factor(latent_k)), piratepal(palette = "pony")[c(4:5)]))
  
  table(data.frame(latent_k.col, Y))
  # Projection de l'ensemble des études.
  plot(X.Y, xlab="PCo1", ylab="PCo2", pch=16, col = Y_col, axes = FALSE)
  legend("bottomright", legend=c("IBD" , "HLT"), fill = unique(Y_col), bty="n", cex=0.8)
  title(paste0("PCoA - ", NAMES[i]))
  
  for (i in unique(latent_k.col)){
    mu <- colMeans(X.Y[which(latent_k.col==i),1:2])
    Sigma <- cov(X.Y[which(latent_k.col==i),1:2]) 
    addEllipse(mu, Sigma, p.interval = 0.8, col = i, lty = 1, lwd=4)
  }
  axis(1)
  axis(2)
  
  # PART 2 ================================================
  plot(x = rep(1, k), y = seq(1,k,by=1),  col = yarrr::piratepal("pony")[4:5], 
       axes = F, xlab = "", cex = 5, pch = 15, ylim = c(0.5,k+.5), ylab = "Alternative clustering")
  for(i in 1:k){text(x = 1, y = i, xpd = NA, c("1", "2","3")[i])}
  for(i in 1:k){text(x = 1, y = i-0.35, xpd = NA,
                     paste(table(Y[which(latent_k==i)])[1],"/", table(Y[which(latent_k==i)])[2]))}
  text(1, k+.5, "Healthy / IBD", xpd = NA)
  text(1, k+.4, "p-value < 0.001", xpd = NA)
  chisq.test(Y, latent_k.col)
  
  comparaison2 = c()
  for(i in 1:k){
    comparaison2 <- cbind(comparaison2, colMeans(res.genus[which(latent_k==i),]))
  }
  composition.genus = comparaison2[w.b,]
  Others = 1 - colSums(comparaison2[w.b,])
  composition.genus = rbind(composition.genus, Others)
  
  ## Barplot
  par(mar=c(4,4,2,2))
  barplot(composition.genus, horiz = T, col = yarrr::piratepal("info2"), 
          main='Gut microbiota composition',  border=NA)
  
  # PART 3 ================================================
  rf <- randomForest(latent_k~., data = X)
  rf$importance
  plot(log(colMeans(X)), log(apply(X,2, var)), 
       cex = rf$importance/2, axes = FALSE, main = "Importance of each species on clustering",
       # xlim=c(-21, 5), ylim=c(-35, 5), 
       xlab = "Species' abundances (log scale)", ylab = "Species' variability (log scale)")
  axis(1)
  axis(2)
  
  dev.off()
}

for(i in 1:length(DF)){
  pdf(file = paste0("figures/figure_1_varImp_", NAMES[i], ".pdf"), width = 5, height = 3.5)
  latent_k <- as.factor(cutree(hclust(DF[[i]], method = "ward.D2"), 2))
  rf <- randomForest(latent_k~., data = X)
  varImpPlot(rf, n=10)
  dev.off()
}

## LEGENDS
pdf(file = paste0("figures/figure_1_legend.pdf"), width = 12.5, height = 3.5)
par(mar=c(0,0,0,0))
plot.new()
legend("center",legend = rownames(composition.genus), title="Genus of microbiome composition",
       fill=yarrr::piratepal("info2"), ncol = 5, bty = "n", cex = 0.7, title.adj = 0)
dev.off()
# legend("bottomright", legend=rep("",3), pch=rep(1,3), pt.cex=c(1,3,5), ncol=4, bty="n",
#        title="Species' importance to define \n latent sample classes", title.cex=0.5)
