BDE_analyze_mt_2 <- function(DATA, MT_Accuracy_set, GENERATION, Heat, AnalyseName, OBJFUNC.ARGS) {
  
  mt_control <- matrix(nrow = 4, ncol = 5)
  colnames(mt_control) <- c('G1', paste0('G',GENERATION), 'all markers', 'first_set', 'final_set')
  rownames(mt_control) <- c('mean_cor_cov', 'sd_cor_cov', 'mean_pheno_cor', 'sd_pheno_cor')
  
  ### Reference correlations
  nTraits <- ncol(DATA$p.probe)
  Y <- DATA$p.probe
  A <- A.mat(as.matrix(DATA$m.probe))
  #
  Rand_effects <- list(
    list(
      K = A, 
      COV = list(
        type = 'UN', 
        df0 = OBJFUNC.ARGS$df0, 
        S0 = diag(nTraits)
      )
    )
  )
  Residuals <- list(
    type = 'DIAG', 
    S0 = rep(OBJFUNC.ARGS$S0, nTraits), 
    df0 = rep(OBJFUNC.ARGS$df1, nTraits)
  )
  capture.output(fm <- MTM(Y = Y,
            K = Rand_effects,
            resCov = Residuals,
            nIter = OBJFUNC.ARGS$nIter,
            burnIn = OBJFUNC.ARGS$burnIn,
            thin = OBJFUNC.ARGS$thin,
            saveAt = OBJFUNC.ARGS$saveAt
  ), file = 'temp/fm')
  cov_ref <- fm$K[[1]]$G
  cor_cov_ref <- cor(cov_ref, DATA$weight_cov)
  cor_cov_accuracy_ref <- mean(sapply(1:ncol(cor_cov_ref), function(i) cor_cov_ref[i,i]))
  prod_accuracy_ref <- cor(fm$YHat[,1], DATA$p.probe[,1])
  
  ### Hit final features correlations
  nTraits <- ncol(DATA$p.probe)
  Y <- DATA$p.probe
  A <- A.mat(as.matrix(DATA$m.probe[,Heat$final_features]))
  #
  Rand_effects <- list(
    list(
      K = A, 
      COV = list(
        type = 'UN', 
        df0 = OBJFUNC.ARGS$df0, 
        S0 = diag(nTraits)
      )
    )
  )
  Residuals <- list(
    type = 'DIAG', 
    S0 = rep(OBJFUNC.ARGS$S0, nTraits), 
    df0 = rep(OBJFUNC.ARGS$df1, nTraits)
  )
  capture.output(fm <- MTM(Y = Y,
            K = Rand_effects,
            resCov = Residuals,
            nIter = OBJFUNC.ARGS$nIter,
            burnIn = OBJFUNC.ARGS$burnIn,
            thin = OBJFUNC.ARGS$thin,
            saveAt = OBJFUNC.ARGS$saveAt
  ), file = 'temp/fm')
  cov_final <- fm$K[[1]]$G
  cor_cov_final <- cor(cov_final, DATA$weight_cov)
  cor_cov_accuracy_final <- mean(sapply(1:ncol(cor_cov_final), function(i) cor_cov_final[i,i]))
  prod_accuracy_final <- cor(fm$YHat[,1], DATA$p.probe[,1])
  
  ### Hit features in G1 correlations
  nTraits <- ncol(DATA$p.probe)
  Y <- DATA$p.probe
  A <- A.mat(as.matrix(DATA$m.probe[,Heat$best_features_in_G1]))
  #
  Rand_effects <- list(
    list(
      K = A, 
      COV = list(
        type = 'UN', 
        df0 = OBJFUNC.ARGS$df0, 
        S0 = diag(nTraits)
      )
    )
  )
  Residuals <- list(
    type = 'DIAG', 
    S0 = rep(OBJFUNC.ARGS$S0, nTraits), 
    df0 = rep(OBJFUNC.ARGS$df1, nTraits)
  )
  capture.output(fm <- MTM(Y = Y,
            K = Rand_effects,
            resCov = Residuals,
            nIter = OBJFUNC.ARGS$nIter,
            burnIn = OBJFUNC.ARGS$burnIn,
            thin = OBJFUNC.ARGS$thin,
            saveAt = OBJFUNC.ARGS$saveAt
  ), file = 'temp/fm')
  cov_first <- fm$K[[1]]$G
  cor_cov_first <- cor(cov_first, DATA$weight_cov)
  cor_cov_accuracy_first <- mean(sapply(1:ncol(cor_cov_first), function(i) cor_cov_first[i,i]))
  prod_accuracy_first <- cor(fm$YHat[,1], DATA$p.probe[,1])
  
  ### Write values in mt_control
  mt_control[1,1] <- mean(MT_Accuracy_set$cor_cov_accuracy_first)
  mt_control[2,1] <- sd(MT_Accuracy_set$cor_cov_accuracy_first)
  mt_control[3,1] <- mean(MT_Accuracy_set$prod_accuracy_first)
  mt_control[4,1] <- sd(MT_Accuracy_set$prod_accuracy_first)
  
  mt_control[1,2] <- mean(MT_Accuracy_set$cor_cov_accuracy_final)
  mt_control[2,2] <- sd(MT_Accuracy_set$cor_cov_accuracy_final)
  mt_control[3,2] <- mean(MT_Accuracy_set$prod_accuracy_final)
  mt_control[4,2] <- sd(MT_Accuracy_set$prod_accuracy_final)
  
  mt_control[1,3] <- cor_cov_accuracy_ref
  mt_control[3,3] <- prod_accuracy_ref
  
  mt_control[1,4] <- cor_cov_accuracy_first
  mt_control[3,4] <- prod_accuracy_first
  mt_control[1,5] <- cor_cov_accuracy_final
  mt_control[3,5] <- prod_accuracy_final
  
  write.csv(mt_control, file = paste0("results/", AnalyseName, "_mt_control.csv"))
  
  ### Join cov_matrices
  upper_cov_matrix <- cbind(DATA$weight_cov, rep(NA, nrow(cov_ref)), cov_ref)
  lower_cov_matrix <- cbind(cov_first, rep(NA, nrow(cov_ref)), cov_final)
  middle_layer <- rep(NA, 2*ncol(cov_ref)+1)
  colnames(upper_cov_matrix) <- c('weight_1', 'weight_2', 'weight_3', '#', 'weight_1', 'weight_2', 'weight_3')
  colnames(lower_cov_matrix) <- c('weight_1', 'weight_2', 'weight_3', '#', 'weight_1', 'weight_2', 'weight_3')
  names(middle_layer) <- c('weight_1', 'weight_2', 'weight_3', 'orig_vs_ref', 'weight_1', 'weight_2', 'weight_3')
  cov_matricies <- rbind(upper_cov_matrix, middle_layer, lower_cov_matrix)
  cov_matricies[4,4] <- 'first_vs_final'
  colnames(cov_matricies) <- c('weight_1', 'weight_2', 'weight_3', 'orig_vs_ref', 'weight_1', 'weight_2', 'weight_3')
  rownames(cov_matricies) <- c('1.weight_1', '2.weight_2', '3.weight_3', '4.###', '5.weight_1', '6.weight_2', '7.weight_3')
  write.csv(cov_matricies, file = paste0("results/", AnalyseName, "_cov_matrices.csv"))
  
  ### Join cor_matrices
  upper_cor_matrix <- cbind(cov2cor(as.matrix(DATA$weight_cov)), rep(NA, nrow(cov_ref)), cov2cor(cov_ref))
  lower_cor_matrix <- cbind(cov2cor(cov_first), rep(NA, nrow(cov_ref)), cov2cor(cov_final))
  middle_layer <- rep(NA, 2*ncol(cov_ref)+1)
  colnames(upper_cor_matrix) <- c('weight_1', 'weight_2', 'weight_3', '#', 'weight_1', 'weight_2', 'weight_3')
  colnames(lower_cor_matrix) <- c('weight_1', 'weight_2', 'weight_3', '#', 'weight_1', 'weight_2', 'weight_3')
  names(middle_layer) <- c('weight_1', 'weight_2', 'weight_3', 'orig_vs_ref', 'weight_1', 'weight_2', 'weight_3')
  cor_matricies <- rbind(upper_cor_matrix, middle_layer, lower_cor_matrix)
  cor_matricies[4,4] <- 'first_vs_final'
  colnames(cor_matricies) <- c('weight_1', 'weight_2', 'weight_3', 'orig_vs_ref', 'weight_1', 'weight_2', 'weight_3')
  rownames(cor_matricies) <- c('1.weight_1', '2.weight_2', '3.weight_3', '4.###', '5.weight_1', '6.weight_2', '7.weight_3')
  write.csv(cor_matricies, file = paste0("results/", AnalyseName, "_cor_matrices.csv"))
}