### Set BDE_parameters
CROSSVAL <- 3
OFFSET <- 0.1
NBASEFEAT <- 400
CFSBEST <- 100 #isNotUsed
NP <- 20
GENERATION <- 15
MUTFACTOR <- 0.3
CR <- 0.5
SEEDRNG <- 12
LMD_CENT <- 0.1
NUMCORES <- ceiling(detectCores()/2)
#
BDE_Parameters <- c(CROSSVAL, OFFSET, NBASEFEAT, CFSBEST, NP, GENERATION, MUTFACTOR, CR, SEEDRNG, LMD_CENT, NUMCORES)
names(BDE_Parameters) <- c('CROSSVAL', 'OFFSET', 'NBASEFEAT', 'CFSBEST', 'NP', 'GENERATION', 'MUTFACTOR', 'CR', 'SEEDRNG', 'LMD_CENT', 'NUMCORES')