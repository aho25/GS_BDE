rm(list=ls())

### Load packages
library(stats)
library(parallel)

### Source scripts
source("BDE_parts/BDE_function.R")
source("gbs_functions/gbs_source.R")
source("BDE_data/data_source.R")
source("BDE_parts/BDE_analyze_synthetic.R")
source("BDE_parts/BDE_analyze_synthetic2.R")

### Set OBJFUNC and validation function parameters  
source("BDE_parameters/BDE_OBJFUNC_parameters.R")

### Set BDE_parameters
source("BDE_parameters/BDE_parameters.R")
OBJFUNC <- gbs_rrblup
AnalyseName <- "rrblup_01"
k <- 0.1 # Choose noise coefficient

### Load data
DATA <- Synthetic_data_ST(k)

set.seed(SEEDRNG)
### Start BDE ###
Accuracy_set <- list()
start_time <- Sys.time()
for (p in 1:5) {
  capture.output(Population <- BDE(DATA$p.probe, DATA$m.probe, CROSSVAL, OFFSET, NBASEFEAT, CFSBEST, NP, GENERATION, MUTFACTOR, CR, SEEDRNG, OBJFUNC, OBJFUNC.ARGS),
                 file = 'temp/Population')
  cat('Population', p, 'done\n')
  Accuracy <- BDE_analyze(DATA, Population)
  Accuracy_set$accuracy_final <- c(Accuracy_set$accuracy_final, Accuracy$accuracy_final)
  Accuracy_set$prod_accuracy_final <- c(Accuracy_set$prod_accuracy_final, Accuracy$prod_accuracy_final)
  Accuracy_set$accuracy_first <- c(Accuracy_set$accuracy_first, Accuracy$accuracy_first)
  Accuracy_set$prod_accuracy_first <- c(Accuracy_set$prod_accuracy_first, Accuracy$prod_accuracy_first)
}
end_time <- Sys.time()
BDE_time <- ceiling(end_time - start_time)

### Analyse ###
BDE_analyze2(DATA, Accuracy_set, GENERATION)
BDE_time

