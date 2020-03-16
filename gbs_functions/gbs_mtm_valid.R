library(MTM)
library(rrBLUP)

### Define Multitrait function for final_features validation
gbs_mtm_valid <- function(PHENO_TRAIN, MARKERS_TRAIN, PHENO_TEST, MARKERS_TEST, VAL.ARGS) {
  nTraits <- ncol(PHENO_TRAIN)
  Y <- rbind(PHENO_TRAIN, matrix(NA, nrow = nrow(PHENO_TEST), ncol = nTraits))
  A <- A.mat(rbind(MARKERS_TRAIN, MARKERS_TEST))
  #
  Rand_effects <- list(
    list(
      K = A, 
      COV = list(
        type = 'UN', 
        df0 = VAL.ARGS$df0, 
        S0 = diag(nTraits)
      )
    )
  )
  Residuals <- list(
    type = 'DIAG', 
    S0 = rep(VAL.ARGS$S0, nTraits), 
    df0 = rep(VAL.ARGS$df1, nTraits)
  )
  fm <- MTM(Y = Y,
            K = Rand_effects,
            resCov = Residuals,
            nIter = VAL.ARGS$nIter,
            burnIn = VAL.ARGS$burnIn,
            thin = VAL.ARGS$thin,
            saveAt = VAL.ARGS$saveAt
  )
  #
  prod_accuracy <- sapply(1:nTraits, function(j) {
    cor(fm$YHat[-c(1:nrow(PHENO_TRAIN)),j], PHENO_TEST[,j])
  })
  return(prod_accuracy)
}