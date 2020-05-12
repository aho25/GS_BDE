source("BDE_parts/BDE_normalize.R")
source("BDE_parts/BDE_optim_parallel.R")
source("BDE_parts/BDE_features_parallel.R")
source("BDE_parts/BDE_IG_16.R")
#source("BDE_parts/BDE_CFS.R")

### Define BDE Algorithm
BDE <- function(PHENO, MARKERS, CROSSVAL, OFFSET, NBASEFEAT, CFSBEST, NP, GENERATION, MUTFACTOR, CR, SEEDRNG, OBJFUNC, OBJFUNC.ARGS, NUMCORES) {
  
  ############### Count feature scores by 3 groups ##################
  Features <- BDE_features(PHENO, MARKERS, OFFSET, NUMCORES)
  feature_names <- colnames(MARKERS)
  
  ### Write scores in vectors
  One_Way_ANOVA <- vector(length = length(feature_names))
  Fisher_Score <- vector(length = length(feature_names))
  for (j in 1:length(feature_names)) {
    One_Way_ANOVA[j] <- Features[[feature_names[j]]]$One_Way_ANOVA$F_value
    Fisher_Score[j] <- Features[[feature_names[j]]]$Fisher_Score
  }
  names(Fisher_Score) <- feature_names
  names(One_Way_ANOVA) <- feature_names
  
  ### Compare One_Way_ANOVA & Fisher_Score results
  One_Way_ANOVA.sorted <- sort(One_Way_ANOVA[which(One_Way_ANOVA != 0)], decreasing = T)
  Fisher_Score.sorted <- sort(Fisher_Score[which(Fisher_Score != 0)], decreasing = T)
  #which(rank(-One_Way_ANOVA.sorted) <= 10)
  #which(rank(-Fisher_Score.sorted) <= 10)
  
  IG.FSelector.sorted_16 <- BDE_IG_16(PHENO, MARKERS)
  
  ### Leave NBASEFEAT best from each filter method
  One_Way_ANOVA.best <- One_Way_ANOVA.sorted[1:NBASEFEAT]
  One_Way_ANOVA.best <- One_Way_ANOVA.best[!is.na(One_Way_ANOVA.best)]
  Fisher_Score.best <- Fisher_Score.sorted[1:NBASEFEAT]
  Fisher_Score.best <- Fisher_Score.best[!is.na(Fisher_Score.best)]
  IG.FSelector_16.best <- IG.FSelector.sorted_16[1:NBASEFEAT]
  IG.FSelector_16.best <- IG.FSelector_16.best[!is.na(IG.FSelector_16.best)]
  
  ### Count overlapping percentage between three filter methods
  One_Way_ANOVA.best_VS_Fisher_Score.best <- length(which(names(One_Way_ANOVA.best) %in% names(Fisher_Score.best)))
  One_Way_ANOVA.best_VS_IG.FSelector_16.best <- length(which(names(One_Way_ANOVA.best) %in% names(IG.FSelector_16.best)))
  Fisher_Score.best_VS_IG.FSelector_16.best <- length(which(names(Fisher_Score.best) %in% names(IG.FSelector_16.best)))
  
  ### Form feature pool
  feature_pool.names <- unique(c(names(One_Way_ANOVA.best), names(Fisher_Score.best), names(IG.FSelector_16.best)))
  
  ### Count Correlation based Feature Selection (CFS) in feature pool
  #feature_pool.names <- BDE_CFS(PHENO, MARKERS, feature_pool.names, CFSBEST)
  
  ### Score normalization
  One_Way_ANOVA.norm <- normalize(One_Way_ANOVA.sorted)
  Fisher_Score.norm <- normalize(Fisher_Score.sorted)
  IG.FSelector_16.norm <- normalize(IG.FSelector.sorted_16)
  
  ### Count probabilities of features
  One_Way_ANOVA.score <- One_Way_ANOVA.norm[feature_pool.names]
  One_Way_ANOVA.score[is.na(One_Way_ANOVA.score)] <- 0
  Fisher_Score.score <- Fisher_Score.norm[feature_pool.names]
  Fisher_Score.score[is.na(Fisher_Score.score)] <- 0
  IG.FSelector_16.score <- IG.FSelector_16.norm[feature_pool.names]
  IG.FSelector_16.score[is.na(IG.FSelector_16.score)] <- 0
  Total_Score <- One_Way_ANOVA.score + Fisher_Score.score
  if (length(IG.FSelector_16.score) != 0) {
    Total_Score <- Total_Score + IG.FSelector_16.score
  }
  #any(is.na(Total_Score))
  #feature_pool.names <- names(Total_Score)
  MARKERS.pool <- as.data.frame(MARKERS[,feature_pool.names])
  p_of_feature <- (max(Total_Score) - Total_Score) / (max(Total_Score) - min(Total_Score))
  D <- length(feature_pool.names) #number of features in feature pool
  for (i in 1:D) {
    p_of_feature[i] <- 1 - min(0.9, max(0.1, p_of_feature[i])) 
  }
  
  ####### BDE ######
  #set.seed(SEEDRNG)
  
  Population <- BDE_optim(PHENO, MARKERS.pool, feature_pool.names, p_of_feature, D, CROSSVAL, NP, GENERATION, MUTFACTOR, CR, SEEDRNG, OBJFUNC, OBJFUNC.ARGS, NUMCORES)
  return(Population)
}
