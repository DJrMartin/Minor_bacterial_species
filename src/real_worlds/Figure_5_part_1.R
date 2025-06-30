rm(list=ls())
library(ape)
library(vegan)
library(yarrr)
library(readxl)
library(cluster)

### data ================
load("data/data.rda")

curated.data <- data$`Wu et al.`
curated.class <- as.factor(data$`Wu et al.`$group) # cela réduit considérablement le jeu de données.
curated.study <- as.factor(data$`Wu et al.`$projectid)

# Colors of the co-variates
col.class <- as.character(factor(curated.class, levels(curated.class), c("#00640090", "#7F00FF90")))
col.study <- as.character(factor(curated.study, levels(curated.study), yarrr::piratepal("pony")[1:length(levels(curated.study))]))
col2rgb("#972C8DFF")[,1]

# Normalisation of the data
tss_data <- data.frame(t(apply(curated.data[,-c(1:4)], 1, function(x) x/sum(x))))
tss_data <- tss_data[,colSums(tss_data)!=0]
clr_data <- compositions::clr(tss_data)

##### Check the composition of gut bacterial ecosystem depending of enterotype =========
taxonomie <- read.table("data/taxonomy_Wu_2024.txt", sep='\t', header=T)
labels=tidyr::separate(data.frame(S=taxonomie$tax_name), 'S', c('K', 'P', 'C', "O", "F", "G", "S"), sep = ";" )
labels = unique(labels)
rownames(labels) = gsub(" ",".",labels$S)

genus = labels[colnames(tss_data),]$G
length(which(is.na(genus)))/length(genus)

tss.genus <- tss_data 
colnames(tss.genus)=genus

tss.G=NULL
for (i in unique(genus)){
  tss.G=cbind(tss.G,apply(as.matrix(tss_data[,which(genus==i)]), 1, sum))
}
colnames(tss.G) = unique(genus)
tss.G <- data.frame(tss.G)

# PCoA
BC <- vegdist(tss_data)
BC_genus <- vegdist(tss.G)
D.clr <- dist(clr_data)
res.pcoa <- pcoa(BC)
res.pcoa.genus <- pcoa(BC_genus)
res.pcoa.D <- pcoa(D.clr)

##### Check the dependence between enterotype and UC =========
layout(matrix(c(1,1,1,7,9,6,6,6,7,7,7,
                2,2,2,3,4,5,5,5,8,8,8,
                2,2,2,3,4,5,5,5,8,8,8,
                2,2,2,3,4,5,5,5,8,8,8), nrow=4, byrow = T))
par(mar=c(2,4,2,1))
# HAC with ward.D method + 
k <- cutree(hclust(BC_genus, method = "ward.D"), 3)
# Panel A
boxplot(res.pcoa.genus$vectors[,1]~curated.class,axes=F, horizontal = T, xlab='', ylab="",
        col=c("#00640090", "#7F00FF90"))
text(x=0.45, y=1.5, paste("log p value =\n",round(log(t.test(res.pcoa.genus$vectors[,1]~curated.class)$p.value)), "(***)"), 
     cex=1, srt=270, xpd=NA)
segments(x0 =0.40, y0=1, y1=2, xpd=NA, lwd=1.5)
# Panel B
par(mar=c(4,4,1,1))
Eigenvalues <- round(res.pcoa.genus$values[1:2,3]*100)
plot(res.pcoa.genus$vectors, xlab=paste0('PCo1 (',Eigenvalues[1],'%)'), ylab=paste0('PCo2 (',Eigenvalues[2],'%)'),
    cex=1, pch=20, col='white')
legend('topleft', legend=c("Healthy", "UC"), fill = c("#00640090", "#7F00FF90"), bty='n', cex=0.9)
# ellipses of enterotype
for (i in unique(k)){
  mu <- colMeans(res.pcoa.genus$vectors[which(k==i),1:2])
  Sigma <- cov(res.pcoa.genus$vectors[which(k==i),1:2]) 
  addEllipse(mu, Sigma, p.interval = 0.8, col = yarrr::piratepal("pony")[i], lty = 1, lwd=4)
}
points(res.pcoa.genus$vectors, col=col.class, pch=20)
# Panel C
par(mar=c(4,2,1,2))
boxplot(res.pcoa.genus$vectors[,2]~curated.class,axes=F, horizontal = F, xlab='', ylab="",
        col=c("#00640090", "#7F00FF90"))
text(x=1.5, y=0.6, paste("log p value =\n",round(log(t.test(res.pcoa.genus$vectors[,2]~curated.class)$p.value)),"(*)"), 
     cex=1, srt=0, xpd=NA)
segments(x0 = 1, x1 = 2, y0=0.55, xpd=NA, lwd=1.5)

comparaison = c()
for(i in unique(k)){
  comparaison <- cbind(comparaison, colMeans(tss.G[which(k==i),]))
}

w <- which(apply(comparaison, 1, var)>0.00017)
length(w)

data = comparaison[w,]
Others = 1 - colSums(comparaison[w,])
data = rbind(data, Others)

# FIGURE 2 : Determine the composition of 
# gut bacterial ecosystem depending of the enterotype
par(mar=c(4,4,2,2))
# dt %>% plot(horiz = T, main='', axes=F)
plot(x=rep(1, length(unique(k))), y = seq(1,length(unique(k)),by=1),  col=yarrr::piratepal("pony")[1:3], 
     axes=F, xlab="", cex=5, pch=15, ylim=c(0.5,3.5), ylab="Enterotypes")
for(i in unique(k)){text(x=1, y=i-0.35, xpd=NA,
                         paste(table(curated.class[which(k==i)])[1],"/", table(curated.class[which(k==i)])[2]))}
text(1, 3.5, "Healthy / UC", xpd=NA)

par(mar=c(4,2,2,2))
barplot(data, horiz = T, col = yarrr::piratepal("info2"), 
        main='Gut microbiota composition',  border=NA)

# Legend of gut bacterial ecosystem
par(mar=c(0,0,0,0))
plot.new()
legend("center",legend=rownames(data), 
       fill=yarrr::piratepal("info2"), ncol=3, bty="n")
plot.new()
legend("right", legend = c("Low", "Medium", "High"), pt.cex=c(0.5,2,4), pch = 1,text.width=0.1,
       title='Species influences on\nenterotype determination', bty="n", ncol=3)

# LAST panel
par(mar=c(4,4,2,2))
rf = randomForest::randomForest(as.factor(k)~., data= tss_data)
plot(log(colMeans(tss_data)), log(apply(tss_data, 2, var)), 
     cex=rf$importance/4, xlab="Species means (log scale)", ylab="Species variance (log scale)")
# legend("bottomright", legend = c("low", "medium", "high"), pt.cex=c(0.5,2,4), pch = 1,
       # title='Species influences on\nenterotype determination', bty="n")

