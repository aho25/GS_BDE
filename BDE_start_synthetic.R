rm(list=ls())

### Load packages
library(stats)
library(parallel)

### Source scripts
source("BDE_parts/BDE_function.R")
source("gbs_functions/gbs_source.R")
source("BDE_data/data_source.R")
source("BDE_parts/BDE_analyze_synthetic.R")

### Set OBJFUNC and validation function parameters  
source("BDE_parameters/BDE_OBJFUNC_parameters.R")

### Set BDE_parameters
source("BDE_parameters/BDE_parameters.R")
OBJFUNC <- gbs_mtm
AnalyseName <- "mtm_synthetic_..."

### Load data
DATA <- Synthetic_data_ST2()

### Start BDE ###
Pops <- list()
start_time <- Sys.time()
for (p in 1:5) {
  capture.output(Population <- BDE(DATA$p.probe, DATA$m.probe, CROSSVAL, OFFSET, NBASEFEAT, CFSBEST, NP, GENERATION, MUTFACTOR, CR, SEEDRNG, OBJFUNC, OBJFUNC.ARGS),
                 file = 'temp/Population')
  cat('Population', p, 'done\n')
  ### Final features
  Pops$final_features[[p]] <- colnames(Population[[paste0('G',GENERATION)]]$X[,which(Population[[paste0('G',GENERATION)]]$X[Population[[paste0('G',GENERATION)]]$x_best,] == 0)])
  ### Best features in first generation
  Pops$features_in_G1[[p]] <- colnames(Population$G1$X[,which(Population$G1$X[Population$G1$x_best,] == 0)])
  ### Best Fitness
  Pops$final_fitness[p] <- Population[[paste0('G',GENERATION)]]$Fitness[Population[[paste0('G',GENERATION)]]$x_best]
}
end_time <- Sys.time()
BDE_time <- ceiling(end_time - start_time)

### Analyse ###
Heat <- BDE_analyze(DATA, Pops, OBJFUNC_Parameters, BDE_Parameters, AnalyseName, BDE_time)
Heat$best_in_G1_heat.MR
Heat$final_heat.MR


