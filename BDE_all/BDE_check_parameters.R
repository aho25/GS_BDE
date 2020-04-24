  ### Set BDE_parameters
  CROSSVAL <- 3
  OFFSET <- 0.1
  NBASEFEAT <- 50
  CFSBEST <- 100 #isNotUsed
  NP <- 12
  GENERATION <- 4
  MUTFACTOR <- 0.3
  CR <- 0.5
  SEEDRNG <- 12
  NUMCORES <- detectCores(logical = F) #isNotUsed
  #
  BDE_Parameters <- c(CROSSVAL, OFFSET, NBASEFEAT, CFSBEST, NP, GENERATION, MUTFACTOR, CR, SEEDRNG)
  names(BDE_Parameters) <- c('CROSSVAL', 'OFFSET', 'NBASEFEAT', 'CFSBEST', 'NP', 'GENERATION', 'MUTFACTOR', 'CR', 'SEEDRNG')
