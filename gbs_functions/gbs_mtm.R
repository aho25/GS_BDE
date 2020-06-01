library(MTM)
library(rrBLUP)

### Define OBJFFUNC for Multitrait Fitness calculation
gbs_mtm <- function(PHENO, MARKERS, OBJFUNC.ARGS, CROSSVAL, SEEDRNG, LMD) {
  nTraits <- ncol(PHENO)
  Markers.nRow <- nrow(MARKERS)
  if (length(sample(Markers.nRow, Markers.nRow%%CROSSVAL)) == 0) {
    splitdata <- split(order(runif(Markers.nRow)), 1:CROSSVAL)
  } else {
    splitdata <- split(order(runif(Markers.nRow))[-sample(Markers.nRow, Markers.nRow%%CROSSVAL)], 1:CROSSVAL)
  }
  fitness <- sapply(1:CROSSVAL, function(i) {
    Y <- PHENO
    Y[splitdata[[i]],1] <- NA
    A <- A.mat(as.matrix(MARKERS))
    #
    Rand_effects <- list(
      list(
        K = A, 
        COV = list(
          type = 'UN', 
          df0 = OBJFUNC.ARGS$df0, 
          S0 = diag(nTraits)
        )
      )
    )
    Residuals <- list(
      type = 'DIAG', 
      S0 = rep(OBJFUNC.ARGS$S0, nTraits), 
      df0 = rep(OBJFUNC.ARGS$df1, nTraits)
    )
    fm <- MTM(Y = Y,
              K = Rand_effects,
              resCov = Residuals,
              nIter = OBJFUNC.ARGS$nIter,
              burnIn = OBJFUNC.ARGS$burnIn,
              thin = OBJFUNC.ARGS$thin,
              saveAt = OBJFUNC.ARGS$saveAt
    )
    #
    prod_accuracy <- sapply(1:nTraits, function(j) {
      cor(fm$YHat[splitdata[[i]],j], PHENO[splitdata[[i]],j])
    })
    return(prod_accuracy[1])
  })
  return(mean(fitness)*(1 - LMD*ncol(MARKERS)))
}