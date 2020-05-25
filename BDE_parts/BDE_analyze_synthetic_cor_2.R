BDE_analyze_cor_2 <- function(DATA, Accuracy_set, GENERATION, Heat, AnalyseName) {
  
  control <- matrix(nrow = 4, ncol = 5)
  colnames(control) <- c('G1', paste0('G',GENERATION), 'rrBLUP', 'first_set', 'final_set')
  rownames(control) <- c('mean_marker_cor', 'sd_marker_cor', 'mean_pheno_cor', 'sd_pheno_cor')
  
  ### Reference correlations# Корелляция с использованием всех фич в регрессии
  prod_model <- mixed.solve(DATA$p.probe, Z=DATA$m.probe)
  accuracy_ref <- cor(DATA$weight_1, prod_model$u) # Control_1 - marker effects correlation
  # 
  prod_g <- prod_model$u
  prod_mu <- prod_model$beta[1]#смещение регрессии
  prod_predicted <- prod_mu + as.matrix(DATA$m.probe) %*% prod_g#predicted pheno
  prod_accuracy_ref <- cor(prod_predicted, DATA$p.probe) # Control_2 - pheno correlation
  
  ### Hit final features correlations# КОР-и на основании фич, по которым строится heatmap
  prod_model <- mixed.solve(DATA$p.probe, Z=DATA$m.probe[,Heat$final_features])
  weight_1_pred <- rep(0, ncol(DATA$m.probe))
  names(weight_1_pred) <- colnames(DATA$m.probe)
  weight_1_pred[Heat$final_features] <- prod_model$u
  accuracy_final <- as.vector(cor(DATA$weight_1, weight_1_pred))
  # 
  prod_g <- prod_model$u
  prod_mu <- prod_model$beta[1]
  prod_predicted <- prod_mu + as.matrix(DATA$m.probe[,Heat$final_features]) %*% prod_g
  prod_accuracy_final <- as.vector(cor(prod_predicted, DATA$p.probe))
  
  ### Hit features in G1 correlations
  prod_model <- mixed.solve(DATA$p.probe, Z=DATA$m.probe[,Heat$best_features_in_G1])
  weight_1_pred <- rep(0, ncol(DATA$m.probe))
  names(weight_1_pred) <- colnames(DATA$m.probe)
  weight_1_pred[Heat$best_features_in_G1] <- prod_model$u
  accuracy_first <- as.vector(cor(DATA$weight_1, weight_1_pred))
  # 
  prod_g <- prod_model$u
  prod_mu <- prod_model$beta[1]
  prod_predicted <- prod_mu + as.matrix(DATA$m.probe[,Heat$best_features_in_G1]) %*% prod_g
  prod_accuracy_first <- as.vector(cor(prod_predicted, DATA$p.probe))

  ### Write values
  control[1,1] <- mean(Accuracy_set$accuracy_first)
  control[2,1] <- sd(Accuracy_set$accuracy_first)
  control[3,1] <- mean(Accuracy_set$prod_accuracy_first)
  control[4,1] <- sd(Accuracy_set$prod_accuracy_first)
  
  control[1,2] <- mean(Accuracy_set$accuracy_final)
  control[2,2] <- sd(Accuracy_set$accuracy_final)
  control[3,2] <- mean(Accuracy_set$prod_accuracy_final)
  control[4,2] <- sd(Accuracy_set$prod_accuracy_final)
  
  control[1,3] <- accuracy_ref 
  control[3,3] <- prod_accuracy_ref
  
  control[1,4] <- accuracy_first
  control[3,4] <- prod_accuracy_first
  control[1,5] <- accuracy_final
  control[3,5] <- prod_accuracy_final
  
  write.csv(control, file = paste0("results/", AnalyseName, "_control.csv"))
}