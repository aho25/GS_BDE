library(BGLR)

### Define Multitrait function for final_features validation
gbs_bglr_mt_valid <- function(PHENO_TRAIN, MARKERS_TRAIN, PHENO_TEST, MARKERS_TEST, VAL.ARGS) {
  nTraits <- ncol(PHENO_TRAIN)
  ETA <- list(list(X = as.matrix(MARKERS_TRAIN), model = 'BayesB'))
  SVD <- svd(PHENO_TRAIN)
  U <- SVD$u
  D <- diag(SVD$d)
  V <- SVD$v
  B <- matrix(nrow = ncol(MARKERS_TRAIN), ncol = ncol(U))
  for(j in 1:nTraits){
    fm <- BGLR(y = U[,j], ETA = ETA, nIter = VAL.ARGS$nIter, burnIn = VAL.ARGS$burnIn, thin = VAL.ARGS$thin,
               saveAt = VAL.ARGS$saveAt, S0 = VAL.ARGS$S0, df0 = VAL.ARGS$df0, verbose = F)
    B[,j] <- fm$ETA[[1]]$b
  }
  BETA <- B %*% D %*% t(V)
  prod_predicted <- as.matrix(MARKERS_TEST) %*% BETA
  prod_accuracy <- sapply(1:nTraits, function(j) {
    cor(prod_predicted[,j], PHENO_TEST[,j])
  })
  return(prod_accuracy)
}