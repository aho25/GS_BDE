library(BGLR)

### Define OBJFFUNC for Multitrait Fitness calculation
gbs_bglr_mt <- function(PHENO, MARKERS, OBJFUNC.ARGS, CROSSVAL, SEEDRNG, LMD) {
  nTraits <- ncol(PHENO)
  Markers.nRow <- nrow(MARKERS)
  if (length(sample(Markers.nRow, Markers.nRow%%CROSSVAL)) == 0) {
    splitdata <- split(order(runif(Markers.nRow)), 1:CROSSVAL)
  } else {
    splitdata <- split(order(runif(Markers.nRow))[-sample(Markers.nRow, Markers.nRow%%CROSSVAL)], 1:CROSSVAL)
  }
  fitness <- sapply(1:CROSSVAL, function(i) {
    pheno_train <- PHENO[-splitdata[[i]],]
    m_train <- as.matrix(MARKERS[-splitdata[[i]],])
    pheno_test <- PHENO[splitdata[[i]],]
    m_test <- as.matrix(MARKERS[splitdata[[i]],])
    #
    ETA <- list(list(X = m_train, model = 'BayesB'))
    SVD <- svd(pheno_train)
    U <- SVD$u
    D <- diag(SVD$d)
    V <- SVD$v
    B <- matrix(nrow = ncol(m_train), ncol = ncol(U))
    for(j in 1:nTraits){
      fm <- BGLR(y = U[,j], ETA = ETA, nIter = OBJFUNC.ARGS$nIter, burnIn = OBJFUNC.ARGS$burnIn, thin = OBJFUNC.ARGS$thin, 
                 saveAt = OBJFUNC.ARGS$saveAt, S0 = OBJFUNC.ARGS$S0, df0 = OBJFUNC.ARGS$df0, verbose = F)
      B[,j] <- fm$ETA[[1]]$b
    }
    BETA <- B %*% D %*% t(V)
    prod_predicted <- m_test %*% BETA
    prod_accuracy <- sapply(1:nTraits, function(j) {
      cor(prod_predicted[,j], pheno_test[,j])
    })
    return(prod_accuracy[1])
  })
  return(mean(fitness)*(1 - LMD*ncol(MARKERS)))
}