BDE_analyze <- function(Population, GENERATION, OBJFUNC_Parameters, VALID_Parameters) {
  ### Best Fitness
  final_fitness.best <- Population[[paste0('G',GENERATION)]]$Fitness[Population[[paste0('G',GENERATION)]]$x_best]
  ### Final features
  final_features.names <- colnames(Population[[paste0('G',GENERATION)]]$X[,which(Population[[paste0('G',GENERATION)]]$X[Population[[paste0('G',GENERATION)]]$x_best,] == 0)])
  final_features.length <- length(final_features.names)
  cat('final_features.length =', final_features.length, '\n')
  ### Best features in first generation
  best_features_in_G1.names <- colnames(Population$G1$X[,which(Population$G1$X[Population$G1$x_best,] == 0)])
  best_features_in_G1.length <- length(best_features_in_G1.names)
  cat('best_features_in_G1.length =', best_features_in_G1.length, '\n')
  ### Number of coincidences of feautures in first and last Generartion
  final_features_VS_best_features_in_G1 <- length(which(final_features.names %in% best_features_in_G1.names))
  cat('final_features_VS_best_features_in_G1 =', final_features_VS_best_features_in_G1, '\n')
  
  write.csv(final_features.names, file = "temp/final_features_BDE.csv")
  write.csv(best_features_in_G1.names, file = "temp/best_features_in_G1_BDE.csv")
  
  ### Plot best Fitness by Generations
  Fitness_best <- vector(length = GENERATION)
  for (G in 1:GENERATION) {
    Fitness_best[G] <- Population[[paste0('G',G)]]$Fitness[Population[[paste0('G',G)]]$x_best]
  }
  plot(Fitness_best, type = 'b', col = 'blue', xlab = 'Generation', ylab = 'Fitness',
       main = 'Improvement of prediction \n by generations')
  
  ### Write BDE_Result ###
  prediction_accuracy_in_G1 <- VALIDFUNC(DATA$p.probe, DATA$m.probe[,best_features_in_G1.names], DATA$p.valid, DATA$m.valid[,best_features_in_G1.names], VAL.ARGS)
  cat('prediction_accuracy_in_G1 =', prediction_accuracy_in_G1, '\n')
  prediction_accuracy_final <- VALIDFUNC(DATA$p.probe, DATA$m.probe[,final_features.names], DATA$p.valid, DATA$m.valid[,final_features.names], VAL.ARGS)
  cat('prediction_accuracy_final =', prediction_accuracy_final, '\n')
  BDE_Result <- as.data.frame(c(final_fitness.best, final_features.length, best_features_in_G1.length, final_features_VS_best_features_in_G1, 
                                prediction_accuracy_in_G1[1], prediction_accuracy_final[1], 'BDE Parameters:', BDE_Parameters, 'OBJ and VAL Parameters:', 
                                OBJFUNC_Parameters, VALID_Parameters), optional = T)
  row.names(BDE_Result) <- c('final_fitness.best', 'final_features.length', 'best_features_in_G1.length', 'final_features_VS_best_features_in_G1',
                             'prediction_accuracy_in_G1', 'prediction_accuracy_final', '###', names(BDE_Parameters), '####', 
                             names(OBJFUNC_Parameters), names(VALID_Parameters))
  write.table(BDE_Result, file = "temp/BDE_Result.txt", col.names = F)
}
