rm(list=ls())

### Load packages
library(stats)
library(parallel)

### Source scripts
source("BDE_parts/BDE_function.R")
source("gbs_functions/gbs_source.R")
source("BDE_data/data_source.R")
source("BDE_parts/BDE_analyze_synthetic_heat.R")
source("BDE_parts/BDE_analyze_synthetic_cor_1.R")
source("BDE_parts/BDE_analyze_synthetic_cor_2.R")
source("BDE_parts/BDE_analyze_synthetic_mt_1.R")
source("BDE_parts/BDE_analyze_synthetic_mt_2.R")

### Set OBJFUNC and validation function parameters  
source("BDE_parameters/BDE_OBJFUNC_parameters.R")

### Set BDE_parameters
source("BDE_parameters/BDE_parameters.R")
OBJFUNC <- gbs_mtm
AnalyseName <- "mtm_1"
k <- 1 # Choose noise coefficient

### Load data
DATA <- Synthetic_data_MT(k)
LMD <- LMD_CENT/ncol(DATA$m.probe)

set.seed(SEEDRNG)
### Start BDE ###
Accuracy_set <- list()
MT_Accuracy_set <- list()
Pops <- list()
start_time <- Sys.time()
for (p in 1:5) {
  capture.output(Population <- BDE(DATA$p.probe, DATA$m.probe, CROSSVAL, OFFSET, NBASEFEAT, CFSBEST, NP, GENERATION, MUTFACTOR, CR, SEEDRNG, OBJFUNC, OBJFUNC.ARGS, LMD),
                 file = 'temp/Population')
  cat('Population', p, 'done\n')
  
  ### Final features
  Pops$final_features[[p]] <- colnames(Population[[paste0('G',GENERATION)]]$X[,which(Population[[paste0('G',GENERATION)]]$X[Population[[paste0('G',GENERATION)]]$x_best,] == 0)])
  ### Best features in first generation
  Pops$features_in_G1[[p]] <- colnames(Population$G1$X[,which(Population$G1$X[Population$G1$x_best,] == 0)])
  ### Best Fitness
  Pops$final_fitness[p] <- Population[[paste0('G',GENERATION)]]$Fitness[Population[[paste0('G',GENERATION)]]$x_best]
  
  ### Collect correlations
  Accuracy <- BDE_analyze_cor_1(DATA, Population)
  Accuracy_set$accuracy_final <- c(Accuracy_set$accuracy_final, Accuracy$accuracy_final)
  Accuracy_set$prod_accuracy_final <- c(Accuracy_set$prod_accuracy_final, Accuracy$prod_accuracy_final)
  Accuracy_set$accuracy_first <- c(Accuracy_set$accuracy_first, Accuracy$accuracy_first)
  Accuracy_set$prod_accuracy_first <- c(Accuracy_set$prod_accuracy_first, Accuracy$prod_accuracy_first)
  
  ### Collect MT correlations
  MT_Accuracy <- BDE_analyze_mt_1(DATA, Population, OBJFUNC.ARGS)
  MT_Accuracy_set$cor_cov_accuracy_final <- c(MT_Accuracy_set$cor_cov_accuracy_final, MT_Accuracy$cor_cov_accuracy_final)
  MT_Accuracy_set$prod_accuracy_final <- c(MT_Accuracy_set$prod_accuracy_final, MT_Accuracy$prod_accuracy_final)
  MT_Accuracy_set$cor_cov_accuracy_first <- c(MT_Accuracy_set$cor_cov_accuracy_first, MT_Accuracy$cor_cov_accuracy_first)
  MT_Accuracy_set$prod_accuracy_first <- c(MT_Accuracy_set$prod_accuracy_first, MT_Accuracy$prod_accuracy_first)
}
end_time <- Sys.time()
BDE_time <- ceiling(end_time - start_time)

### Analyse ###
Heat <- BDE_analyze_heat(DATA, Pops, OBJFUNC_Parameters, BDE_Parameters, AnalyseName, BDE_time)
BDE_analyze_cor_2(DATA, Accuracy_set, GENERATION, Heat, AnalyseName)
BDE_analyze_mt_2(DATA, MT_Accuracy_set, GENERATION, Heat, AnalyseName, OBJFUNC.ARGS)

Heat$best_in_G1_heat.MR
Heat$final_heat.MR

