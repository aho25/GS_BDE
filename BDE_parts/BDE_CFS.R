library(lsr)

BDE_CFS <- function(PHENO, MARKERS, feature_pool.names, CFSBEST) {
  k <- length(feature_pool.names)
  MARKERS.pool <- as.data.frame(MARKERS[,feature_pool.names])
  CFS_score <- vector(length = k)
  names(CFS_score) <- feature_pool.names
  class_corr <- vector(length = k)
  inter_corr <- matrix(nrow = k, ncol = k)
  inter_corr.av <- vector(length = k) #mean value of each feature from inter_corr
  for (i in 1:k) {
    class_corr[i] <- cor(MARKERS.pool[,i], PHENO, method = "spearman")
  }
  class_corr.av <- mean(class_corr)
  for (i in 1:k) {
    for (j in 1:k) {
      tbl <- table(MARKERS.pool[,i], MARKERS_tot.sc[,j])
      inter_corr[i,j] <- cramersV(tbl)
    }
    inter_corr.av[i] <- mean(inter_corr[i,])
    CFS_score[i] <- k*class_corr.av / sqrt(k + k*(k-1)*inter_corr.av[i])
  }
  CFS_score.best <- cutoff.k(as.data.frame(CFS_score), CFSBEST)
  CFS_score.best <- CFS_score.best[!is.na(CFS_score.best)]
  feature_pool.names <- names(CFS_score.best)
  
  return(feature_pool.names)
}
