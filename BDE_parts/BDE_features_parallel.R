BDE_features <- function(PHENO, MARKERS, OFFSET, NUMCORES) {
  
  ### Create a list of Features
  feature_names <- colnames(MARKERS)
  Features <- list()
  
  ### Create feature_score_dimnames list to name feature_score matrices
  feature_score_dimnames <- list()
  feature_score_dimnames$score <- c('average_value', 'variance', 'standard_deviation', 'number_of_samples', 'Sum', 'sum_of_squares')
  feature_score_dimnames$groups <- c('negative', 'null', 'positive')
  
  Features <- mcmapply(function(j) {
    #	cat(j,feature_names[j], '\n')
    Features[[feature_names[j]]]$yield <- list(feature_score_dimnames$groups) #Form 3 groups of feature_yield data
    Features[[feature_names[j]]]$score <- matrix(data = NA, nrow = 6, ncol = 3, dimnames = feature_score_dimnames) #Form feature_score matrix
    ### Write feature_yield data by groups
    for (i in 1:nrow(MARKERS)) {
      if (MARKERS[i,j] < -OFFSET) {
        Features[[feature_names[j]]]$yield$negative <- c(Features[[feature_names[j]]]$yield$negative, PHENO[i,1])
        MARKERS[i,j] <- -1
      } else if (MARKERS[i,j] > OFFSET) {
        Features[[feature_names[j]]]$yield$positive <- c(Features[[feature_names[j]]]$yield$positive, PHENO[i,1])
        MARKERS[i,j] <- 1
      } else {
        Features[[feature_names[j]]]$yield$null <- c(Features[[feature_names[j]]]$yield$null, PHENO[i,1])
        MARKERS[i,j] <- 0
      }
    }
    ### Count and write scores ('average_value', 'variance', 'standard_deviation', 'number_of_samples', 'Sum', 'sum_of_squares') in feature_score matrix:
    ### for 'negative' group
    if (is.null(Features[[feature_names[j]]]$yield$negative) == F) {
      Features[[feature_names[j]]]$score[1,1] <- mean(Features[[feature_names[j]]]$yield$negative)
      Features[[feature_names[j]]]$score[2,1] <- var(Features[[feature_names[j]]]$yield$negative)
      Features[[feature_names[j]]]$score[3,1] <- sd(Features[[feature_names[j]]]$yield$negative)
      Features[[feature_names[j]]]$score[4,1] <- length(Features[[feature_names[j]]]$yield$negative)
      Features[[feature_names[j]]]$score[5,1] <- sum(Features[[feature_names[j]]]$yield$negative)
      Features[[feature_names[j]]]$score[6,1] <- sum((Features[[feature_names[j]]]$yield$negative)^2)
    }
    ### for 'null' group
    if (is.null(Features[[feature_names[j]]]$yield$null) == F) {
      Features[[feature_names[j]]]$score[1,2] <- mean(Features[[feature_names[j]]]$yield$null)
      Features[[feature_names[j]]]$score[2,2] <- var(Features[[feature_names[j]]]$yield$null)
      Features[[feature_names[j]]]$score[3,2] <- sd(Features[[feature_names[j]]]$yield$null)
      Features[[feature_names[j]]]$score[4,2] <- length(Features[[feature_names[j]]]$yield$null)
      Features[[feature_names[j]]]$score[5,2] <- sum(Features[[feature_names[j]]]$yield$null)
      Features[[feature_names[j]]]$score[6,2] <- sum((Features[[feature_names[j]]]$yield$null)^2)
    }
    ### for 'positive' group
    if (is.null(Features[[feature_names[j]]]$yield$positive) == F) {
      Features[[feature_names[j]]]$score[1,3] <- mean(Features[[feature_names[j]]]$yield$positive)
      Features[[feature_names[j]]]$score[2,3] <- var(Features[[feature_names[j]]]$yield$positive)
      Features[[feature_names[j]]]$score[3,3] <- sd(Features[[feature_names[j]]]$yield$positive)
      Features[[feature_names[j]]]$score[4,3] <- length(Features[[feature_names[j]]]$yield$positive)
      Features[[feature_names[j]]]$score[5,3] <- sum(Features[[feature_names[j]]]$yield$positive)
      Features[[feature_names[j]]]$score[6,3] <- sum((Features[[feature_names[j]]]$yield$positive)^2)
    }
    
    ### Count average value of yield by feature
    Features[[feature_names[j]]]$average_yield <- mean(c(Features[[feature_names[j]]]$score[1,1], Features[[feature_names[j]]]$score[1,2],
                                                         Features[[feature_names[j]]]$score[1,3]))
    
    ### Count One_Way_ANOVA scores
    Features[[feature_names[j]]]$One_Way_ANOVA$N <- sum(Features[[feature_names[j]]]$score[4,1], Features[[feature_names[j]]]$score[4,2],
                                                        Features[[feature_names[j]]]$score[4,3]) #Total number of samples used in feature_score counts
    Features[[feature_names[j]]]$One_Way_ANOVA$SSt <- sum(c((Features[[feature_names[j]]]$yield$negative - Features[[feature_names[j]]]$average_yield)^2,
                                                            (Features[[feature_names[j]]]$yield$null - Features[[feature_names[j]]]$average_yield)^2,
                                                            (Features[[feature_names[j]]]$yield$positive - Features[[feature_names[j]]]$average_yield)^2)) #Total sum of squares (SSt = SSbw + SSwg)
    Features[[feature_names[j]]]$One_Way_ANOVA$SSwg <- sum(c((Features[[feature_names[j]]]$yield$negative - Features[[feature_names[j]]]$score[1,1])^2,
                                                             (Features[[feature_names[j]]]$yield$null - Features[[feature_names[j]]]$score[1,2])^2,
                                                             (Features[[feature_names[j]]]$yield$positive - Features[[feature_names[j]]]$score[1,3])^2)) #Sum of squares within groups
    Features[[feature_names[j]]]$One_Way_ANOVA$SSbg <- sum(c(Features[[feature_names[j]]]$score[4,1]*(Features[[feature_names[j]]]$score[1,1] - Features[[feature_names[j]]]$average_yield)^2,
                                                             Features[[feature_names[j]]]$score[4,2]*(Features[[feature_names[j]]]$score[1,2] - Features[[feature_names[j]]]$average_yield)^2,
                                                             Features[[feature_names[j]]]$score[4,3]*(Features[[feature_names[j]]]$score[1,3] - Features[[feature_names[j]]]$average_yield)^2)) #Sum of squares between groups
    Features[[feature_names[j]]]$One_Way_ANOVA$MSwg <- Features[[feature_names[j]]]$One_Way_ANOVA$SSwg / (Features[[feature_names[j]]]$One_Way_ANOVA$N - 3) #Mean of squares within groups
    Features[[feature_names[j]]]$One_Way_ANOVA$MSbg <- Features[[feature_names[j]]]$One_Way_ANOVA$SSbg / 2 #Mean of squares between groups
    Features[[feature_names[j]]]$One_Way_ANOVA$F_value <- Features[[feature_names[j]]]$One_Way_ANOVA$MSbg / Features[[feature_names[j]]]$One_Way_ANOVA$MSwg #Result score of the test
    
    ### Count Fisher_Score
    Features[[feature_names[j]]]$Fisher_Score <- sum(c(Features[[feature_names[j]]]$score[4,1]*(Features[[feature_names[j]]]$score[1,1] - Features[[feature_names[j]]]$average_yield)^2,
                                                       Features[[feature_names[j]]]$score[4,2]*(Features[[feature_names[j]]]$score[1,2] - Features[[feature_names[j]]]$average_yield)^2,
                                                       Features[[feature_names[j]]]$score[4,3]*(Features[[feature_names[j]]]$score[1,3] - Features[[feature_names[j]]]$average_yield)^2)
    ) / sum(c(Features[[feature_names[j]]]$score[4,1]*Features[[feature_names[j]]]$score[2,1]^2,
              Features[[feature_names[j]]]$score[4,2]*Features[[feature_names[j]]]$score[2,2]^2,
              Features[[feature_names[j]]]$score[4,3]*Features[[feature_names[j]]]$score[2,3]^2))
    return(Features)
  }, 1:length(feature_names), mc.cores = NUMCORES)
  
  return(Features)
}
