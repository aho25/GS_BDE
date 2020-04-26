library(BMTME)

### Define OBJFFUNC for Multitrait Fitness calculation
gbs_bme <- function(PHENO, MARKERS, OBJFUNC.ARGS, CROSSVAL, SEEDRNG, NUMCORES) {
  PHENO <- as.data.frame(PHENO)
  MARKERS <- as.matrix(MARKERS)
  A <- tcrossprod(MARKERS)/ncol(MARKERS)
  LG <- cholesky(A)
  ZG <- model.matrix(~0 + row.names(PHENO))
  Z.G <- ZG %*% LG
  Y <- as.matrix(PHENO)
  #
  pheno <- data.frame(GID = row.names(PHENO), Response = PHENO$yield)
  CrossV <- CV.KFold(pheno, DataSetID = 'GID', K = CROSSVAL, set_seed = SEEDRNG)
  #
  pm <- BME(Y = Y, Z1 = Z.G, nIter = OBJFUNC.ARGS$nIter, burnIn = OBJFUNC.ARGS$burnIn, thin = OBJFUNC.ARGS$thin, bs = OBJFUNC.ARGS$bs, parallelCores = NUMCORES, testingSet = CrossV)
  #
  accuracy <- summary(pm)
  return(accuracy$Pearson[which(as.vector(accuracy$Trait) == 'yield')])
}