library(BGLR)

### Define OBJFFUNC for Multitrait Fitness calculation
gbs_bglr_st <- function(PHENO, MARKERS, OBJFUNC.ARGS, CROSSVAL, SEEDRNG) {
  Markers.nRow <- nrow(MARKERS)
  if (length(sample(Markers.nRow, Markers.nRow%%CROSSVAL)) == 0) {
    splitdata <- split(order(runif(Markers.nRow)), 1:CROSSVAL)
  } else {
    splitdata <- split(order(runif(Markers.nRow))[-sample(Markers.nRow, Markers.nRow%%CROSSVAL)], 1:CROSSVAL)
  }
  fitness <- sapply(1:CROSSVAL, function(i) {
    pheno_train <- PHENO[-splitdata[[i]],1]
    m_train <- as.matrix(MARKERS[-splitdata[[i]],])
    pheno_test <- PHENO[splitdata[[i]],1]
    m_test <- as.matrix(MARKERS[splitdata[[i]],])
    #
    ETA <- list(list(X = m_train, model = 'BayesB'))
    fm <- BGLR(y = pheno_train, ETA = ETA, nIter = OBJFUNC.ARGS$nIter, burnIn = OBJFUNC.ARGS$burnIn, thin = OBJFUNC.ARGS$thin, 
               saveAt = OBJFUNC.ARGS$saveAt, S0 = OBJFUNC.ARGS$S0, df0 = OBJFUNC.ARGS$df0, verbose = F)
    prod_predicted <- m_test %*% fm$ETA[[1]]$b
    prod_accuracy <- cor(prod_predicted, pheno_test)
    return(prod_accuracy)
  })
  return(mean(fitness))
}