library("heatmaply")
BDE_analyze <- function(DATA, Population, GENERATION, OBJFUNC_Parameters, BDE_Parameters, AnalyseName, BDE_time) {
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
  ### Plot best Fitness by Generations
  Fitness_best <- vector(length = GENERATION)
  for (G in 1:GENERATION) {
    Fitness_best[G] <- Population[[paste0('G',G)]]$Fitness[Population[[paste0('G',G)]]$x_best]
  }
  png(paste0("results/", AnalyseName, "_Fithess_by_Generation.png"))
  plot(Fitness_best, type = 'b', col = 'blue', xlab = 'Generation', ylab = 'Fitness', main = 'Improvement of prediction \n by generations')
  dev.off()
  
  ### Final features rate
  correct_final_features.length <- length(which(final_features.names %in% DATA$major_snp_1)) #количество правильно предсказанных фич
  correct_final_features.rate <- correct_final_features.length / length(final_features.names) #количество правильных к количеству предсказанных
  cat('\n correct_final_features.rate =', correct_final_features.rate, '\n')
  error_final_features.length <- length(final_features.names) - correct_final_features.length #ошибочнo предсказанных
  error_final_features.rate <- error_final_features.length / length(final_features.names) 
  cat('error_final_features.rate =', error_final_features.rate, '\n')
  caught_final_features.length <- correct_final_features.length
  caught_final_features.rate <- caught_final_features.length / length(DATA$major_snp_1)
  cat('caught_final_features.rate =', caught_final_features.rate, '\n')
  lost_final_features.length <- length(DATA$major_snp_1) - caught_final_features.length
  lost_final_features.rate <- lost_final_features.length / length(DATA$major_snp_1)
  cat('lost_final_features.rate =', lost_final_features.rate, '\n')
  final_features.MR <- matrix(data = c(correct_final_features.length, 
                                       lost_final_features.length, 
                                       error_final_features.length,
                                       (ncol(DATA$m.probe) - correct_final_features.length - lost_final_features.length - error_final_features.length)), 
                              ncol = 2)
  colnames(final_features.MR) <- c("Major", "Rest")
  rownames(final_features.MR) <- c("Major", "Rest")
  final_heat.MR <- heatmaply(final_features.MR, dendrogram = "none", hide_colorbar = T, draw_cellnote = T, colors = rainbow(100, alpha = 0.75, start = 0, end = 0.15, rev = T),
                             cellnote_size = 22, grid_color = "white", column_text_angle = 0, xlab = 'Actual SNP', ylab = 'Predicted SNP',
                             main = 'Final features    .', fontsize_row = 12, fontsize_col = 12, subplot_widths = 0.4, subplot_heights = 1, 
                             cellnote_textposition = "middle center", cellnote_color = "black", colorbar_xpos = 3, 
                             file = paste0("results/", AnalyseName, "_final_heat.MR.html"))
  
  ### Best features in G1 rate
  correct_features_in_G1.length <- length(which(final_features.names %in% DATA$major_snp_1))
  correct_features_in_G1.rate <- correct_features_in_G1.length / length(final_features.names)
  cat('\n correct_features_in_G1.rate =', correct_features_in_G1.rate, '\n')
  error_features_in_G1.length <- length(final_features.names) - correct_features_in_G1.length
  error_features_in_G1.rate <- error_features_in_G1.length / length(final_features.names)
  cat('error_features_in_G1.rate =', error_features_in_G1.rate, '\n')
  caught_features_in_G1.length <- correct_features_in_G1.length
  caught_features_in_G1.rate <- caught_features_in_G1.length / length(DATA$major_snp_1)
  cat('caught_features_in_G1.rate =', caught_features_in_G1.rate, '\n')
  lost_features_in_G1.length <- length(DATA$major_snp_1) - caught_features_in_G1.length
  lost_features_in_G1.rate <- lost_features_in_G1.length / length(DATA$major_snp_1)
  cat('lost_features_in_G1.rate =', lost_features_in_G1.rate, '\n')
  features_in_G1.MR <- matrix(data = c(correct_features_in_G1.length, lost_features_in_G1.length, error_features_in_G1.length,
                                       (ncol(DATA$m.probe) - correct_features_in_G1.length - lost_features_in_G1.length - error_features_in_G1.length)), ncol = 2)
  colnames(features_in_G1.MR) <- c("Major", "Rest")
  rownames(features_in_G1.MR) <- c("Major", "Rest")
  best_in_G1_heat.MR <- heatmaply(features_in_G1.MR, dendrogram = "none", hide_colorbar = T, draw_cellnote = T, colors = rainbow(100, alpha = 0.75, start = 0, end = 0.15, rev = T),
                                  cellnote_size = 22, grid_color = "white", column_text_angle = 0, xlab = 'Actual SNP', ylab = 'Predicted SNP',
                                  main = 'Features in G1     .', fontsize_row = 12, fontsize_col = 12, subplot_widths = 0.4, subplot_heights = 1, 
                                  cellnote_textposition = "middle center", cellnote_color = "black", colorbar_xpos = 3, 
                                  file = paste0("results/", AnalyseName, "_best_in_G1_heat.MR.html"))
  
  ### Write BDE_Result ###
  BDE_time_print <- capture.output(BDE_time)
  BDE_Result <- as.data.frame(c(final_fitness.best, final_features.length, correct_final_features.rate, error_final_features.rate,
                                caught_final_features.rate, lost_final_features.rate, best_features_in_G1.length, 
                                correct_features_in_G1.rate, error_features_in_G1.rate, caught_features_in_G1.rate, 
                                lost_features_in_G1.rate, final_features_VS_best_features_in_G1, 
                                BDE_time_print, 'BDE Parameters:', BDE_Parameters, 'OBJ_FUNC Parameters:', 
                                OBJFUNC_Parameters), optional = T)
  row.names(BDE_Result) <- c('final_fitness.best', 'final_features.length', 'correct_final_features.rate', 'error_final_features.rate',
                             'caught_final_features.rate', 'lost_final_features.rate', 'best_features_in_G1.length', 
                             'correct_features_in_G1.rate', 'error_features_in_G1.rate', 'caught_features_in_G1.rate', 
                             'lost_features_in_G1.rate', 'final_features_VS_best_features_in_G1',
                             'BDE_time', '###', names(BDE_Parameters), '####', names(OBJFUNC_Parameters))
  write.table(BDE_Result, file = paste0("results/", AnalyseName, "_BDE_Result.txt"), col.names = F)
  return(list(final_heat.MR=final_heat.MR, best_in_G1_heat.MR=best_in_G1_heat.MR))
}
