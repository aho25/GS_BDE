rm(list=ls())
#mtm + mtm_valid with Soy_MT
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
AnalyseName <- "mtm_synthetic_data..."

### Load data
DATA <- Synthetic_data_MT(OFFSET)

### Start BDE ###
start_time <- Sys.time()
system.time(
  Population <- BDE(DATA$p.probe, DATA$m.probe, CROSSVAL, OFFSET, NBASEFEAT, CFSBEST, NP, GENERATION, MUTFACTOR, CR, SEEDRNG, OBJFUNC, OBJFUNC.ARGS)
)
end_time <- Sys.time()
BDE_time <- ceiling(end_time - start_time)

### Analyse ###
Heat <- BDE_analyze(DATA, Population, GENERATION, OBJFUNC_Parameters, BDE_Parameters, AnalyseName, BDE_time)
Heat$best_in_G1_heat.MR
Heat$final_heat.MR
  

