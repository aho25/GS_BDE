library(FSelector)

BDE_IG_16 <- function(PHENO, MARKERS) {
  ### Split 'yield' data by 16 classes & form a dataframe to count IG
  PHENO.max <- max(PHENO)
  PHENO.min <- min(PHENO) - 1
  PHENO.seq_16 <- (PHENO.max - PHENO.min) / 16
  #hist(PHENO, breaks = seq(PHENO.min, PHENO.max, PHENO.seq_16), col = 'orange')
  PHENO.classes_16 <- cut(PHENO, breaks = seq(PHENO.min, PHENO.max, PHENO.seq_16), labels = letters[1:16])
  MARKERS.classes_16 <- cbind(MARKERS, as.vector(PHENO.classes_16))
  colnames(MARKERS.classes_16) <- c(colnames(MARKERS), 'classes')
  MARKERS.classes_16 <- as.data.frame(MARKERS.classes_16)
  
  ### Count IG & check results
  IG.FSelector_16 <- information.gain(classes ~ ., data=MARKERS.classes_16)
  IG.FSelector.vector_16 <- as.vector(IG.FSelector_16[,1])
  names(IG.FSelector.vector_16) <- rownames(IG.FSelector_16)
  IG.FSelector.sorted_16 <- sort(IG.FSelector.vector_16[which(IG.FSelector.vector_16 != 0)], decreasing = T)
  #which(rank(-IG.FSelector.sorted_16) <= 10)
  
  return(IG.FSelector.sorted_16)
}
