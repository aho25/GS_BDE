BDE_analyze_mt_1 <- function(DATA, Population, OBJFUNC.ARGS) {
  
  ### Final features
  final_features.names <- colnames(Population[[paste0('G',GENERATION)]]$X[,which(Population[[paste0('G',GENERATION)]]$X[Population[[paste0('G',GENERATION)]]$x_best,] == 0)])
  ### Best features in first generation
  features_in_G1.names <- colnames(Population$G1$X[,which(Population$G1$X[Population$G1$x_best,] == 0)])
  
  ### Correlations_in_final_generation
  nTraits <- ncol(DATA$p.probe)
  Y <- DATA$p.probe
  A <- A.mat(as.matrix(DATA$m.probe[,final_features.names]))
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
  
  ### Correlations_in_first_generation
  nTraits <- ncol(DATA$p.probe)
  Y <- DATA$p.probe
  A <- A.mat(as.matrix(DATA$m.probe[,features_in_G1.names]))
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
  
  return(list(cor_cov_accuracy_final=cor_cov_accuracy_final, prod_accuracy_final=prod_accuracy_final,
              cor_cov_accuracy_first=cor_cov_accuracy_first, prod_accuracy_first=prod_accuracy_first))
}
  