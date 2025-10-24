rm(list=ls())

### Packages.
packages <- readLines("requirements.txt")
invisible(lapply(packages, function(pkg){library(pkg, character.only = TRUE)}))

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
  
### Test for the number of cross-validation runs
set.seed(2)

MARTIN = PRJEB1220 = list()
for(i in seq(5, 20, 5)){
  # MARTIN et al.
  res = out_of_bag_prediction(X, Y, transformed = FALSE, cv = i, p = 0.8)
  MARTIN[[paste("nb_", i)]]  = cbind(res$AUC.X, res$AUC.R)
  
  # H. B. Nielsen, et al.
  res = out_of_bag_prediction(X_2, Y_2, transformed = FALSE, cv = i, p = 0.8)
  PRJEB1220[[paste("nb_", i)]] = cbind(res$AUC.X, res$AUC.R)
  # Print the results
  print(i)
}

color = rep(NA, 6)
color[seq(1,5,by=2)] = LaCroixColoR::lacroix_palette(type = "paired")[seq(2,6,by=2)]
color[seq(2,6,by=2)] = LaCroixColoR::lacroix_palette(type = "paired")[seq(1,5,by=2)]

pdf(file = "figures/figure_S5A.pdf", width = 8.5, height = 4)
#### PLOT concerning the number of cross-validation.
layout(matrix(c(1,2), nrow=1))
# PLOT A (PRJEB1220)
df = lapply(PRJEB1220, function(x){apply(x, 2, mean)})
sd_res <- lapply(PRJEB1220, function(x){apply(x, 2, sd)})

plot(x = 0.9:3.9, unlist(df)[seq(1, 7, 2)], type = 'b', ylim = c(0.5, 1), xlim = c(0.8, 4.2), 
     ylab = 'AUROC', xlab = "Number of cross-validations", axes = F, col = color[1], pch = 20)
points(x = 1.1:4.1, unlist(df)[seq(2, 8, 2)], type = 'b', ylim = c(0.5, 1), col = color[2], pch = 20)
segments(x0 = 0.9:3.9, y0 = unlist(df)[seq(1, 7, 2)]-unlist(sd_res)[seq(1, 7, 2)], y1 = unlist(df)[seq(1, 7, 2)]+unlist(sd_res)[seq(1, 7, 2)], col = color[1])
segments(x0 = 1.1:4.1, y0 = unlist(df)[seq(2, 8, 2)]-unlist(sd_res)[seq(2, 8, 2)], y1 = unlist(df)[seq(2, 8, 2)]+unlist(sd_res)[seq(2, 8, 2)], col = color[2])
axis(1, at = 1:4, labels = c(5, 10, 15, 20))
axis(2)
text(x = 1, y = 0.5, "PRJEB1220", pos = 4, cex = 0.8)
# PLOT B (MARTIN)
df = lapply(MARTIN, function(x){apply(x, 2, mean)})
sd_res <- lapply(MARTIN, function(x){apply(x, 2, sd)})

plot(x = 0.9:3.9, unlist(df)[seq(1, 7, 2)], type = 'b', ylim = c(0.5, 1), xlim = c(0.8, 4.2), 
     ylab = 'AUROC', xlab = "Number of cross-validations", axes = F, col = color[1], pch = 20)
points(x = 1.1:4.1, unlist(df)[seq(2, 8, 2)], type = 'b', ylim = c(0.5, 1), col = color[2], pch = 20)
segments(x0 = 0.9:3.9, y0 = unlist(df)[seq(1, 7, 2)]-unlist(sd_res)[seq(1, 7, 2)], y1 = unlist(df)[seq(1, 7, 2)]+unlist(sd_res)[seq(1, 7, 2)], col = color[1])
segments(x0 = 1.1:4.1, y0 = unlist(df)[seq(2, 8, 2)]-unlist(sd_res)[seq(2, 8, 2)], y1 = unlist(df)[seq(2, 8, 2)]+unlist(sd_res)[seq(2, 8, 2)], col = color[2])
axis(1, at = 1:4, labels = c(5, 10, 15, 20))
axis(2)
text(x = 1, y = 0.5, "Martin et al. (2025)", pos = 4, cex = 0.8)

legend("top", legend = c("rel. ab.", "rel. ab. transformed"), fill = color[1:2], cex = 0.8, bty = "n")
dev.off()

#### Tests concerning the split between train-tests
set.seed(2)
MARTIN = WU = c()
for(split in seq(0.6, 0.9, by = 0.1)){
  # MARTIN et al.
  res = out_of_bag_prediction(X, Y, transformed = FALSE, cv = 20, p = split)
  MARTIN[[paste("nb_", split)]]  = cbind(res$AUC.X, res$AUC.R)
  
  # H. B. Nielsen, et al.
  res = out_of_bag_prediction(X_2, Y_2, transformed = FALSE, cv = 20, p = split)
  PRJEB1220[[paste("nb_", split)]] = cbind(res$AUC.X, res$AUC.R)
  
  # Print the number of the for()
  print(split)
}

pdf(file = "figures/figure_S5B.pdf", width = 8.5, height = 4)
#### PLOT concerning the number of cross-validation.
layout(matrix(c(1,2), nrow=1))
# PLOT A (PRJEB1220)
df = lapply(PRJEB1220, function(x){apply(x, 2, mean)})
sd_res <- lapply(PRJEB1220, function(x){apply(x, 2, sd)})

plot(x = 0.9:3.9, unlist(df)[seq(1, 7, 2)], type = 'b', ylim = c(0.5, 1), xlim = c(0.8, 4.2), 
     ylab = 'AUROC', xlab = "Proportion of samples in train set", axes = F, col = color[1], pch = 20)
points(x = 1.1:4.1, unlist(df)[seq(2, 8, 2)], type = 'b', ylim = c(0.5, 1), col = color[2], pch = 20)
segments(x0 = 0.9:3.9, y0 = unlist(df)[seq(1, 7, 2)]-unlist(sd_res)[seq(1, 7, 2)], y1 = unlist(df)[seq(1, 7, 2)]+unlist(sd_res)[seq(1, 7, 2)], col = color[1])
segments(x0 = 1.1:4.1, y0 = unlist(df)[seq(2, 8, 2)]-unlist(sd_res)[seq(2, 8, 2)], y1 = unlist(df)[seq(2, 8, 2)]+unlist(sd_res)[seq(2, 8, 2)], col = color[2])
axis(1, at = 1:4, labels = seq(0.6, 0.9, by = 0.1))
axis(2)
text(x = 1, y = 0.5, "PRJEB1220", pos = 4, cex = 0.8)
# PLOT B (MARTIN)
df = lapply(MARTIN, function(x){apply(x, 2, mean)})
sd_res <- lapply(MARTIN, function(x){apply(x, 2, sd)})

plot(x = 0.9:3.9, unlist(df)[seq(1, 7, 2)], type = 'b', ylim = c(0.5, 1), xlim = c(0.8, 4.2), 
     ylab = 'AUROC', xlab = "Proportion of samples in train set", axes = F, col = color[1], pch = 20)
points(x = 1.1:4.1, unlist(df)[seq(2, 8, 2)], type = 'b', ylim = c(0.5, 1), col = color[2], pch = 20)
segments(x0 = 0.9:3.9, y0 = unlist(df)[seq(1, 7, 2)]-unlist(sd_res)[seq(1, 7, 2)], y1 = unlist(df)[seq(1, 7, 2)]+unlist(sd_res)[seq(1, 7, 2)], col = color[1])
segments(x0 = 1.1:4.1, y0 = unlist(df)[seq(2, 8, 2)]-unlist(sd_res)[seq(2, 8, 2)], y1 = unlist(df)[seq(2, 8, 2)]+unlist(sd_res)[seq(2, 8, 2)], col = color[2])
axis(1, at = 1:4, labels = seq(0.6, 0.9, by = 0.1))
axis(2)
text(x = 1, y = 0.5, "Martin et al. (2025)", pos = 4, cex = 0.8)

legend("top", legend = c("rel. ab.", "rel. ab. transformed"), fill = color[1:2], cex = 0.8, bty = "n")
dev.off()

#### Tests concerning the hyper parameters.
set.seed(2)
MARTIN = PRJEB1220 = list()
for(maxnodes in c(5, 25, 50)){
  res.PRJEB1220 = res.martin = list()
  for(mtry in c(5, 25, 50)){
    # MARTIN et al.
    res = out_of_bag_prediction(X, Y, transformed = FALSE, cv = 10,
                                p = 0.7, maxnodes = maxnodes, mtry = mtry)
    res.martin[[paste("mtry", mtry)]] = cbind(res$AUC.X, res$AUC.R)
    # WU et al.
    res = out_of_bag_prediction(X_2, Y_2, transformed = FALSE, cv = 10, 
                                p = 0.7, maxnodes = maxnodes, mtry = mtry)
    res.PRJEB1220[[paste("mtry", mtry)]] = cbind(res$AUC.X, res$AUC.R)
    
    print(mtry)
  }
  
  MARTIN[[paste('maxnodes =', maxnodes)]] = res.martin
  PRJEB1220[[paste('maxnodes =', maxnodes)]] = res.PRJEB1220
  print(maxnodes)
}

pdf(file = "figures/figure_S5C.pdf", width = 8.5, height = 4)
layout(matrix(c(1,2), nrow=1))
# PLOT A (PRJEB1220)
df_5 = lapply(PRJEB1220$`maxnodes = 5`, function(x){apply(x, 2, mean)})
sd_res_5 <- lapply(PRJEB1220$`maxnodes = 5`, function(x){apply(x, 2, sd)})
df_25 = lapply(PRJEB1220$`maxnodes = 25`, function(x){apply(x, 2, mean)})
sd_res_25 <- lapply(PRJEB1220$`maxnodes = 25`, function(x){apply(x, 2, sd)})
df_50 = lapply(PRJEB1220$`maxnodes = 50`, function(x){apply(x, 2, mean)})
sd_res_50 <- lapply(PRJEB1220$`maxnodes = 50`, function(x){apply(x, 2, sd)})

# maxnodes = 5
plot(x = 0.7:3.7, unlist(df_5)[seq(1, 7, 2)], type = 'b', ylim = c(0.5, 1), xlim = c(0.8, 4.2), 
     ylab = 'AUROC', xlab = "mtry", axes = F, col = color[1], pch = 20)
points(x = 1.1:4.1, unlist(df_5)[seq(2, 8, 2)], type = 'b', ylim = c(0.5, 1), col = color[2], pch = 20)
# maxnodes = 25
points(x = 0.8:3.8, unlist(df_25)[seq(1, 7, 2)], type = 'b', ylim = c(0.5, 1), col = color[1], pch = 20, lty = 2)
points(x = 1.2:4.2, unlist(df_25)[seq(2, 8, 2)], type = 'b', ylim = c(0.5, 1), col = color[2], pch = 20, lty = 2)
# maxnodes = 50
points(x = 0.9:3.9, unlist(df_50)[seq(1, 7, 2)], type = 'b', ylim = c(0.5, 1), col = color[1], pch = 20, lty = 3)
points(x = 1.3:4.3, unlist(df_50)[seq(2, 8, 2)], type = 'b', ylim = c(0.5, 1), col = color[2], pch = 20, lty = 3)
segments(x0 = 0.7:3.7, y0 = unlist(df_5)[seq(1, 7, 2)]-unlist(sd_res_5)[seq(1, 7, 2)], y1 = unlist(df_5)[seq(1, 7, 2)]+unlist(sd_res_5)[seq(1, 7, 2)], col = color[1])
segments(x0 = 1.1:4.1, y0 = unlist(df_5)[seq(2, 8, 2)]-unlist(sd_res_5)[seq(2, 8, 2)], y1 = unlist(df_5)[seq(2, 8, 2)]+unlist(sd_res_5)[seq(2, 8, 2)], col = color[2])
segments(x0 = 0.8:3.8, y0 = unlist(df_25)[seq(1, 7, 2)]-unlist(sd_res_25)[seq(1, 7, 2)], y1 = unlist(df_25)[seq(1, 7, 2)]+unlist(sd_res_25)[seq(1, 7, 2)], col = color[1])
segments(x0 = 1.2:4.2, y0 = unlist(df_25)[seq(2, 8, 2)]-unlist(sd_res_25)[seq(2, 8, 2)], y1 = unlist(df_25)[seq(2, 8, 2)]+unlist(sd_res_25)[seq(2, 8, 2)], col = color[2])
segments(x0 = 0.9:3.9, y0 = unlist(df_50)[seq(1, 7, 2)]-unlist(sd_res_50)[seq(1, 7, 2)], y1 = unlist(df_50)[seq(1, 7, 2)]+unlist(sd_res_50)[seq(1, 7, 2)], col = color[1])
segments(x0 = 1.3:4.3, y0 = unlist(df_50)[seq(2, 8, 2)]-unlist(sd_res_50)[seq(2, 8, 2)], y1 = unlist(df_50)[seq(2, 8, 2)]+unlist(sd_res_50)[seq(2, 8, 2)], col = color[2])
axis(1, at = 1:3, labels = c(5, 25, 50))
axis(2)
text(x = 1, y = 0.5, "PRJEB1220", pos = 4, cex = 0.8)

# PLOT B (MARTIN)
df_5 = lapply(MARTIN$`maxnodes = 5`, function(x){apply(x, 2, mean)})
sd_res_5 <- lapply(MARTIN$`maxnodes = 5`, function(x){apply(x, 2, sd)})
df_25 = lapply(MARTIN$`maxnodes = 25`, function(x){apply(x, 2, mean)})
sd_res_25 <- lapply(MARTIN$`maxnodes = 25`, function(x){apply(x, 2, sd)})
df_50 = lapply(MARTIN$`maxnodes = 50`, function(x){apply(x, 2, mean)})
sd_res_50 <- lapply(MARTIN$`maxnodes = 50`, function(x){apply(x, 2, sd)})

# maxnodes = 5
plot(x = 0.7:3.7, unlist(df_5)[seq(1, 7, 2)], type = 'b', ylim = c(0.5, 1), xlim = c(0.8, 4.2), 
     ylab = 'AUROC', xlab = "mtry", axes = F, col = color[1], pch = 20)
points(x = 1.1:4.1, unlist(df_5)[seq(2, 8, 2)], type = 'b', ylim = c(0.5, 1), col = color[2], pch = 20)
# maxnodes = 25
points(x = 0.8:3.8, unlist(df_25)[seq(1, 7, 2)], type = 'b', ylim = c(0.5, 1), col = color[1], pch = 20, lty = 2)
points(x = 1.2:4.2, unlist(df_25)[seq(2, 8, 2)], type = 'b', ylim = c(0.5, 1), col = color[2], pch = 20, lty = 2)
# maxnodes = 50
points(x = 0.9:3.9, unlist(df_50)[seq(1, 7, 2)], type = 'b', ylim = c(0.5, 1), col = color[1], pch = 20, lty = 3)
points(x = 1.3:4.3, unlist(df_50)[seq(2, 8, 2)], type = 'b', ylim = c(0.5, 1), col = color[2], pch = 20, lty = 3)
segments(x0 = 0.7:3.7, y0 = unlist(df_5)[seq(1, 7, 2)]-unlist(sd_res_5)[seq(1, 7, 2)], y1 = unlist(df_5)[seq(1, 7, 2)]+unlist(sd_res_5)[seq(1, 7, 2)], col = color[1])
segments(x0 = 1.1:4.1, y0 = unlist(df_5)[seq(2, 8, 2)]-unlist(sd_res_5)[seq(2, 8, 2)], y1 = unlist(df_5)[seq(2, 8, 2)]+unlist(sd_res_5)[seq(2, 8, 2)], col = color[2])
segments(x0 = 0.8:3.8, y0 = unlist(df_25)[seq(1, 7, 2)]-unlist(sd_res_25)[seq(1, 7, 2)], y1 = unlist(df_25)[seq(1, 7, 2)]+unlist(sd_res_25)[seq(1, 7, 2)], col = color[1])
segments(x0 = 1.2:4.2, y0 = unlist(df_25)[seq(2, 8, 2)]-unlist(sd_res_25)[seq(2, 8, 2)], y1 = unlist(df_25)[seq(2, 8, 2)]+unlist(sd_res_25)[seq(2, 8, 2)], col = color[2])
segments(x0 = 0.9:3.9, y0 = unlist(df_50)[seq(1, 7, 2)]-unlist(sd_res_50)[seq(1, 7, 2)], y1 = unlist(df_50)[seq(1, 7, 2)]+unlist(sd_res_50)[seq(1, 7, 2)], col = color[1])
segments(x0 = 1.3:4.3, y0 = unlist(df_50)[seq(2, 8, 2)]-unlist(sd_res_50)[seq(2, 8, 2)], y1 = unlist(df_50)[seq(2, 8, 2)]+unlist(sd_res_50)[seq(2, 8, 2)], col = color[2])
axis(1, at = 1:3, labels = c(5, 25, 50))
axis(2)
text(x = 1, y = 0.5, "Martin et al. (2025)", pos = 4, cex = 0.8)

legend("top", legend = c(5, 25, 50), lty = c(1:3), cex = 0.8, bty = "n")
dev.off()


