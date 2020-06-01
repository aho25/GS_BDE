rm(list=ls())

### Load packages
library(stats)
library(parallel)

### Source scripts
source("BDE_parts/BDE_function_parallel.R")
source("gbs_functions/gbs_source.R")
source("BDE_data/data_source.R")
source("BDE_parts/BDE_analyze.R")

### Set OBJFUNC and validation function parameters  
source("BDE_parameters/BDE_OBJFUNC_parameters.R")

### Set BDE_parameters
source("BDE_parameters/BDE_parameters_parallel.R")
OBJFUNC <- gbs_rrblup
VALIDFUNC <- gbs_rrblup_valid
AnalyseName <- "rrblup_synthetic_par2...."

### Load data
DATA <- Soy_init_ST()

### Start BDE ###
start_time <- Sys.time()
system.time(
  Population <- BDE(DATA$p.probe, DATA$m.probe, CROSSVAL, OFFSET, NBASEFEAT, CFSBEST, NP, GENERATION, MUTFACTOR, CR, SEEDRNG, OBJFUNC, OBJFUNC.ARGS, NUMCORES)
)
end_time <- Sys.time()
BDE_time <- ceiling(end_time - start_time)

### Analyse ###
BDE_analyze(DATA, Population, GENERATION, OBJFUNC_Parameters, VALID_Parameters, BDE_Parameters, AnalyseName, BDE_time)
