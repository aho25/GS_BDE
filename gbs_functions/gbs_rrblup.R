library(rrBLUP)

### Define OBJFFUNC for Fitness calculation
gbs_rrblup <- function(PHENO, MARKERS, OBJFUNC.ARGS, CROSSVAL, SEEDRNG) {
  Markers.nRow <- nrow(MARKERS)
  if (length(sample(Markers.nRow, Markers.nRow%%CROSSVAL)) == 0) {
    splitdata <- split(order(runif(Markers.nRow)), 1:CROSSVAL)
  } else {
    splitdata <- split(order(runif(Markers.nRow))[-sample(Markers.nRow, Markers.nRow%%CROSSVAL)], 1:CROSSVAL)
  }
  fitness <- mcmapply(function(i) {
    pheno_train <- PHENO[-splitdata[[i]]]
    m_train <- MARKERS[-splitdata[[i]],]
    pheno_test <- PHENO[splitdata[[i]]]
    m_test <- MARKERS[splitdata[[i]],]
    prod_model <- mixed.solve(pheno_train, Z = m_train, K = NULL, SE = FALSE, return.Hinv = FALSE)
    prod_g <- prod_model$u
    prod_mu <- prod_model$beta[1]
    #
    prod_predicted <- prod_mu + as.matrix(m_test) %*% prod_g
    prod_accuracy <- cor(prod_predicted, pheno_test)
    return(prod_accuracy)
  }, 1:CROSSVAL)
  return(mean(fitness))
}
