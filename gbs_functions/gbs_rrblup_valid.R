library(rrBLUP)

### Define function for final_features validation
gbs_rrblup_valid <- function(PHENO_TRAIN, MARKERS_TRAIN, PHENO_TEST, MARKERS_TEST, VAL.ARGS) {
  prod_model <- mixed.solve(PHENO_TRAIN[,1], Z = MARKERS_TRAIN, K = NULL, SE = FALSE, return.Hinv = FALSE)
  prod_g <- prod_model$u
  prod_mu <- prod_model$beta[1]
  #
  prod_predicted <- prod_mu + as.matrix(MARKERS_TEST) %*% prod_g
  prod_accuracy <- cor.test(prod_predicted, PHENO_TEST[,1])
  return(prod_accuracy$estimate)
}