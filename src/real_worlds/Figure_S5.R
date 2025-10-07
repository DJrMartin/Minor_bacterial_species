rm(list=ls())
######## TESTS
### data ================
load("data/data.rda")
source("functions/out_of_bag.R")

#==================================================================================
### MARTIN ET AL. (2025) ###
## X
X <- as.data.frame(t(apply(data$`Martin et al.`$count, 2, function(x) x/sum(x))))
X_2 <- as.data.frame(t(apply(data$PRJEB1220$count, 2, function(x) x/sum(x))))
# rowSums(X)
## Y
FO <- data$`Martin et al.`$data$Lipides_SV1/9
Y <- as.factor(FO>0.4)
Y_2 <- as.factor(data$PRJEB1220$data$disease)
  
### Test for the number of cross validation 
set.seed(2)
MARTIN = WU = list()
for(i in seq(5, 20, 5)){
  res.wu = res.martin = NULL
  for(j in 1:10){
    # MARTIN et al.
    res = out_of_bag_prediction(X, Y, transformed = FALSE, cv = j, p = 0.8)
    res.martin = rbind(res.martin, c(mean(res$AUC.X), mean(res$AUC.R)))
    # MARTIN et al.
    res = out_of_bag_prediction(X_2, Y_2, transformed = FALSE, cv = j, p = 0.8)
    res.wu = rbind(res.wu, c(mean(res$AUC.X), mean(res$AUC.R)))
  }
  MARTIN[[paste('nb_',i)]] = res.martin
  WU[[paste('nb_',i)]] = res.wu
}

color = rep(NA, 6)
color[seq(1,5,by=2)] = LaCroixColoR::lacroix_palette(type = "paired")[seq(2,6,by=2)]
color[seq(2,6,by=2)] = LaCroixColoR::lacroix_palette(type = "paired")[seq(1,5,by=2)]

layout(matrix(c(1,2), nrow=1))
# PLOT A (MARTIN)
df = lapply(MARTIN, function(x){apply(x, 2, mean)})
sd_res <- lapply(MARTIN, function(x){apply(x, 2, sd)})

plot(x = 0.9:3.9, unlist(df)[seq(1, 7, 2)], type = 'b', ylim = c(0.5, 0.9), xlim = c(0.8, 4.2), 
     ylab = 'AUROC', xlab = "Number of cross-validations", axes = F, col = color[1], pch = 20)
points(x = 1.1:4.1, unlist(df)[seq(2, 8, 2)], type = 'b', ylim = c(0.5, 0.9), col = color[2], pch = 20)
segments(x0 = 0.9:3.9, y0 = unlist(df)[seq(1, 7, 2)]-unlist(sd_res)[seq(1, 7, 2)], y1 = unlist(df)[seq(1, 7, 2)]+unlist(sd_res)[seq(1, 7, 2)], col = color[1])
segments(x0 = 1.1:4.1, y0 = unlist(df)[seq(2, 8, 2)]-unlist(sd_res)[seq(2, 8, 2)], y1 = unlist(df)[seq(2, 8, 2)]+unlist(sd_res)[seq(2, 8, 2)], col = color[2])
axis(1, at = 1:4, labels = c(5, 10, 15, 20))
axis(2)
text(x = 1, y = 0.5, "Martin et al. (2025)", pos = 4, cex = 0.8)

# PLOT B (WU)
df = lapply(WU, function(x){apply(x, 2, mean)})
sd_res <- lapply(WU, function(x){apply(x, 2, sd)})

plot(x = 0.9:3.9, unlist(df)[seq(1, 7, 2)], type = 'b', ylim = c(0.5, 0.9), xlim = c(0.8, 4.2), 
     ylab = 'AUROC', xlab = "Number of cross-validations", axes = F, col = color[1], pch = 20)
points(x = 1.1:4.1, unlist(df)[seq(2, 8, 2)], type = 'b', ylim = c(0.5, 0.9), col = color[2], pch = 20)
segments(x0 = 0.9:3.9, y0 = unlist(df)[seq(1, 7, 2)]-unlist(sd_res)[seq(1, 7, 2)], y1 = unlist(df)[seq(1, 7, 2)]+unlist(sd_res)[seq(1, 7, 2)], col = color[1])
segments(x0 = 1.1:4.1, y0 = unlist(df)[seq(2, 8, 2)]-unlist(sd_res)[seq(2, 8, 2)], y1 = unlist(df)[seq(2, 8, 2)]+unlist(sd_res)[seq(2, 8, 2)], col = color[2])
axis(1, at = 1:4, labels = c(5, 10, 15, 20))
axis(2)
text(x = 1, y = 0.5, "Wu et al. (2024)", pos = 4, cex = 0.8)



rm(sd)
### Test for the number of cross validation 
