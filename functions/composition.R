
genus.compo <- function(x, tree){
  labels = tidyr::separate(data.frame(S = tree$tip.label), 'S', c("CAG",'K', 'P', 'C', "O", "F", "G", "S"), sep = "__" )
  tree$tip.label = labels$S
  rownames(labels)=labels$S
  filt_tree <- ape::keep.tip(tree, intersect(tree$tip.label, colnames(X)))
  
  NEW.OTU <- x
  GENUS=labels[colnames(x),]$G
  colnames(NEW.OTU)=labels[rownames(x),]$G
  
  encapsuler=NULL
  for (i in unique(GENUS)){
    encapsuler=cbind(encapsuler,apply(as.matrix(NEW.OTU[,which(GENUS==i)]), 1, sum))
  }
  colnames(encapsuler)=unique(GENUS)
  X <- data.frame(encapsuler)
  
  return(list(X, filt_tree))
}


