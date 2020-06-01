BDE_analyze_cor_1 <- function(DATA, Population) {
  
  ### Final features
  final_features.names <- colnames(Population[[paste0('G',GENERATION)]]$X[,which(Population[[paste0('G',GENERATION)]]$X[Population[[paste0('G',GENERATION)]]$x_best,] == 0)])
  ### Best features in first generation
  features_in_G1.names <- colnames(Population$G1$X[,which(Population$G1$X[Population$G1$x_best,] == 0)])
  
  ### Correlations_final_generation
  prod_model <- mixed.solve(DATA$p.probe[,1], Z=DATA$m.probe[,final_features.names])
  weight_1_pred <- rep(0, ncol(DATA$m.probe))
  names(weight_1_pred) <- colnames(DATA$m.probe)
  weight_1_pred[final_features.names] <- prod_model$u
  accuracy_final <- as.vector(cor(DATA$weight_1, weight_1_pred)) # Control_1 - marker effects correlation
  # 
  prod_g <- prod_model$u
  prod_mu <- prod_model$beta[1]
  prod_predicted <- prod_mu + as.matrix(DATA$m.probe[,final_features.names]) %*% prod_g
  prod_accuracy_final <- as.vector(cor(prod_predicted, DATA$p.probe[,1])) # Control_2 - pheno correlation
  
  ### Correlations_first_generation
  prod_model <- mixed.solve(DATA$p.probe[,1], Z=DATA$m.probe[,features_in_G1.names])
  weight_1_pred <- rep(0, ncol(DATA$m.probe))
  names(weight_1_pred) <- colnames(DATA$m.probe)
  weight_1_pred[features_in_G1.names] <- prod_model$u
  accuracy_first <- as.vector(cor(DATA$weight_1, weight_1_pred)) # Control_1 - marker effects correlation
  # 
  prod_g <- prod_model$u
  prod_mu <- prod_model$beta[1]
  prod_predicted <- prod_mu + as.matrix(DATA$m.probe[,features_in_G1.names]) %*% prod_g
  prod_accuracy_first <- as.vector(cor(prod_predicted, DATA$p.probe[,1])) # Control_2 - pheno correlation
  
  return(list(accuracy_final=accuracy_final, prod_accuracy_final=prod_accuracy_final, accuracy_first=accuracy_first, prod_accuracy_first=prod_accuracy_first))
}
