library(BGLR)

### Define Multitrait function for final_features validation
gbs_bglr_st_valid <- function(PHENO_TRAIN, MARKERS_TRAIN, PHENO_TEST, MARKERS_TEST, VAL.ARGS) {
  ETA <- list(list(X = as.matrix(MARKERS_TRAIN), model = 'BayesB'))
  fm <- BGLR(y = PHENO_TRAIN[,1], ETA = ETA, nIter = VAL.ARGS$nIter, burnIn = VAL.ARGS$burnIn, thin = VAL.ARGS$thin,
             saveAt = VAL.ARGS$saveAt, S0 = VAL.ARGS$S0, df0 = VAL.ARGS$df0, verbose = F)
  prod_predicted <- as.matrix(MARKERS_TEST) %*% fm$ETA[[1]]$b
  prod_accuracy <- cor(prod_predicted, PHENO_TEST[,1])
  return(prod_accuracy)
}