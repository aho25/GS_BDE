rm(list=ls())
#bglr_mt + bglr_mt_valid with Soy_MT
### Load packages
library(stats)
library(parallel)

### Source scripts
source("BDE_parts/BDE_function.R")
source("gbs_functions/gbs_source.R")
source("BDE_data/data_source.R")
source("BDE_all/BDE_analyze.R")

### Load data
DATA <- Soy_init_MT()

### Set OBJFUNC and validation function parameters
source("BDE_all/BDE_OBJFUNC_parameters.R")

### Set BDE_parameters
source("BDE_all/BDE_parameters.R")
OBJFUNC <- gbs_bglr_mt
VALIDFUNC <- gbs_bglr_mt_valid
AnalyseName <- "bglr_mt_valid_with_Soy_MT"

### Start BDE ###
start_time <- Sys.time()
system.time(
  Population <- BDE(DATA$p.probe, DATA$m.probe, CROSSVAL, OFFSET, NBASEFEAT, CFSBEST, NP, GENERATION, MUTFACTOR, CR, SEEDRNG, OBJFUNC, OBJFUNC.ARGS)
)
end_time <- Sys.time()
BDE_time <- ceiling(end_time - start_time)

### Analyse ###
BDE_analyze(Population, GENERATION, OBJFUNC_Parameters, VALID_Parameters, AnalyseName, BDE_time)


