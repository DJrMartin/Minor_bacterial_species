rm(list=ls())
library(ape)
library(vegan)
library(yarrr)
library(readxl)
library(randomForest)

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

### Out-of-sampling
n <- nrow(tss_data)
p <- ncol(tss_data)
k=3
n.CV <- 20

roc.r = roc.CoDa = c()
imp.r = imp.coda = c()
inter.r = inter.coda =c()
set.seed(1)
for(i in 1:n.CV){
  intraining <- caret::createDataPartition(curated.class, p=0.8, list = F)
  TRAIN <- tss_data[intraining,]
  TEST <- tss_data[-intraining,]
  
  # Z_1 estimation
  bc <- vegdist(tss.G)
  k_est <- cutree(hclust(bc, method="ward.D"), k)
  
  # Residuals train
  TRAIN.R <- data.frame((model<-lm(as.matrix(tss_data[intraining,])~0+k_est[intraining]))$residuals)
  # TRAIN.R <- data.frame((model<-lm(as.matrix(TRAIN)~k_est))$residuals)
  # TEST.R <- data.frame((model<-lm(as.matrix(TEST)~k_est.test))$residuals)
  TEST.R = TEST
  for(j in unique(k_est)){
    TEST.R[which(k_est[-intraining]==j),]=scale(TEST[which(k_est[-intraining]==j),],
                                                center = model$coefficients[1,], scale=F)
  }
  rf.r <- randomForest(curated.class[intraining]~., TRAIN.R, maxnodes = 5, mtry=20, ntree=2000)
  rf.CoDa <- randomForest(curated.class[intraining]~., TRAIN, maxnodes = 5, mtry=20, ntree=2000)
  
  imp.r <- cbind(imp.r, rf.r$importance)
  imp.coda <- cbind(imp.coda, rf.CoDa$importance)
  
  roc.r = c(roc.r , pROC::auc(curated.class[-intraining],predict(rf.r, TEST.R, type='prob')[,1]))
  roc.CoDa = c(roc.CoDa, pROC::auc(curated.class[-intraining],predict(rf.CoDa, TEST, type='prob')[,1]))

  ### Interpretability
  w <- which(rf.r$importance!=0)
  inter.r <- c(inter.r, summary(lm(log(colMeans(tss_data))[w]~log(rf.r$importance[w])))$r.squared)
  w <- which(rf.CoDa$importance!=0)
  inter.coda <- c(inter.coda, summary(lm(log(colMeans(tss_data))[w]~log(rf.CoDa$importance[w])))$r.squared)
  
  print(paste("Done", i))
}

color = LaCroixColoR::lacroix_palette(type = "paired")[c(2,1)]

par(mar=c(4,4,3,2))
plot(1, 1, xlim=c(0,1), ylim=c(0,1), col="white",
     xlab="Sensibility", ylab="Specificity",main="Confident Interval of the predictive performance")
legend("bottomleft",legend=c("Relative abundance data", "Guided-transformed data"), title=paste("95% of Confident Interval"),
       fill=color, bty="n")
abline(1,-1)

bc= vegan::vegdist(tss_data)
R <- data.frame((model<-lm(as.matrix(tss_data)~0+cutree(hclust(bc, method="ward.D"), 3)))$residuals)

diff(table(tss_data[,grep("Roseburia.fae", colnames(R))]>0, as.numeric(curated.class)-1)[,2])

n= sample(1:1300, 1000)
roc_obj <- pROC::roc(curated.class[-n], predict(randomForest(curated.class[n]~., tss_data[n,]),tss_data[-n,], type="prob")[,2], ci=TRUE)  
roc_obj.r <- pROC::roc(curated.class[-n], predict(randomForest(curated.class[n]~., R[n,]),R[-n,], type="prob")[,2], ci=TRUE)  

ci_se <- pROC::ci.se(roc_obj, boot.n=2000)
plot(ci_se, type="shape", col="#FC688250", add=T)  # Affiche l'IC sous forme d'ombre
lines(roc_obj, col='#FC6882')

ci_se <- pROC::ci.se(roc_obj.r, boot.n=2000)
plot(ci_se, type="shape", col="#C70E7B50", add=T)  # Affiche l'IC sous forme d'ombre
lines(roc_obj.r, col='#C70E7B')

layout(matrix(c(1,2), nrow=1, byrow = T))
## PANEL B
par(mar=c(4,5,3,2))
boxplot(roc.CoDa, roc.r, axes=F, ylim=c(0.8, 1.03), main = "Predictive performances",
        ylab='AUC-ROC (cross-validations with\n80% of samples in the train set)', xlab="", col=color)
abline(h=0.5, lty=2)
box()
axis(2)
text(x=1.5, 1.02,paste("p value =",round(t.test(roc.r, roc.CoDa)$p.value, 2)))
segments(x0=1,x1= 2,y0=1, lwd=1.5)
legend("bottom",legend=c("Relative abundance data", "Guided-transformed data"),
       fill=color, bty="n", cex=0.8)

## PANEL C
boxplot(inter.coda, inter.r, axes=F, main = "Bacterial signals" ,ylim = c(0, 0.6),
        xlab="", col=color, ylab="Strength of the relation between\nimportance and abundance of species")

t.test(inter.r, inter.coda)
box()
axis(2)
text(x=1.5, 0.55,"***",cex = 1.5)
segments(x0=1,x1= 2,y0=0.53, lwd=1.5)


res.ML <- list(species_importance_r = imp.r, species_importance_coda = imp.coda)
# save(res.ML, file='~/Dropbox/Res_species_importance.RDATA')
load('~/Dropbox/Res_species_importance.RDATA')
network.r <- data.frame(tss_data[,names(sort(apply(res.ML$species_importance_r, 1, sum), decreasing=T)[1:20])], UC=as.numeric(curated.class)-1)
network.coda <- data.frame(tss_data[,names(sort(apply(res.ML$species_importance_coda, 1, sum), decreasing=T)[1:20])], Y=as.numeric(curated.class)-1)

plot(ape::pcoa(dist(R[,names(sort(apply(res.ML$species_importance_r, 1, sum), decreasing=T)[1:20])]))$vectors, )

PV = FROM = TO = Direction = NULL
for(i in 1:dim(network.r)[2]){
  for(j in 2:dim(network.r)[2]){
    if(i!=j){
      PV =c(PV,chisq.test(network.r[,i]>0, network.r[,j]>0)$p.value)
      FROM = c(FROM, colnames(network.r)[i])
      TO = c(TO, colnames(network.r)[j])
      Direction = c(Direction, diff(table(network.r[,i]>0, network.r[,j]>0)[,2]))
    }
  }
}
bacteria <- data.frame(from = FROM,
                       to = TO,
                       Direction = Direction,
                       weight = abs(log(PV)))
library(dplyr)
library(igraph)
bacteria_sorted <- bacteria %>% 
  rowwise() %>%  # Appliquer les opérations ligne par ligne
  mutate(pair = paste(sort(c(from, to)), collapse = "-")) %>%  # Créer une clé unique pour chaque paire
  ungroup() %>%  # Sortir du mode "rowwise"
  distinct(pair, .keep_all = TRUE) %>%  # Garder une seule occurrence par paire
  select(-pair)  # Supprimer la colonne temporaire 'pair'

net <- bacteria_sorted[bacteria_sorted$weight>quantile(bacteria_sorted$weight,0.8),]
net$weight = net$Direction
net <- net %>% select(-Direction)

# Créer le graphe
g <- graph_from_data_frame(net, directed = FALSE)
E(g)$color <- ifelse(E(g)$weight < 0, "firebrick", "cornflowerblue")
E(g)$weight = abs(E(g)$weight/100)

taxonomie <- read.table("taxonomy_ncbi.txt", sep='\t', header=T)
labels=tidyr::separate(data.frame(S=taxonomie$tax_name), 'S', c('K', 'P', 'C', "O", "F", "G", "S"), sep = ";" )
labels = unique(labels)
rownames(labels) = gsub(" ",".",labels$S)
genus = labels[V(g)$name,]$G
col.genus <- match(genus,names(w))
col.genus[is.na(col.genus)]=14
V(g)$color <- yarrr::piratepal("info2")[col.genus]
V(g)$color[16] = "#972C8DFF"
# Définir la mise en page
# layout <- layout_with_fr(g)  # Algorithme de Fruchterman-Reingold pour une bonne répartition
# Visualiser le graphe
plot(g, 
     layout = layout, 
     vertex.size = 13, 
     vertex.label.cex = 0.7, 
     vertex.label.color = "black",
     edge.width = E(g)$weight, 
     edge.color = E(g)$color, 
     main = "Networks of the 20 most important features derived from RF")

plot.new()
legend("center", legend=c("Positive co-occurence", "Negative co-occurence"), col=c("cornflowerblue", "firebrick"), bty="n", lty=1, lwd=5)
plot(res.pcoa.genus$vectors, xlab=paste0('PCo1 (',Eigenvalues[1],'%)'), ylab=paste0('PCo2 (',Eigenvalues[2],'%)'),
     cex=1, pch=20, col=as.factor(rowSums(tss_data[,V(g)$name[which(layout[,2]>1.8)]]>0)>5))

