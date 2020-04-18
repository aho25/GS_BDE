rm(list=ls())

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

### Set BDE_parameters
#source("BDE_all/BDE_parameters.R")
source("BDE_all/BDE_check_parameters.R")
OBJFUNC <- gbs_bglr_mt
VALIDFUNC <- gbs_bglr_mt_valid

### Set OBJFUNC and validation function parameters
source("BDE_all/BDE_OBJFUNC_parameters.R")

### Start BDE ###
system.time(
  Population <- BDE(DATA$p.probe, DATA$m.probe, CROSSVAL, OFFSET, NBASEFEAT, CFSBEST, NP, GENERATION, MUTFACTOR, CR, SEEDRNG, OBJFUNC, OBJFUNC.ARGS)
)

### Analyse ###
BDE_analyze(Population, GENERATION, OBJFUNC_Parameters, VALID_Parameters)


