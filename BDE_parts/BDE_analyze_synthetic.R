library(heatmaply)
BDE_analyze <- function(DATA, Pops, OBJFUNC_Parameters, BDE_Parameters, AnalyseName, BDE_time) {
  ### Find frequent features
  final_features_unique <- vector()
  features_in_G1_unique <- vector()
  for (p in 1:length(Pops$final_features)) {
    final_features_unique <- c(final_features_unique, Pops$final_features[[p]])
    features_in_G1_unique <- c(features_in_G1_unique, Pops$features_in_G1[[p]])
  }
  final_features_unique <- unique(final_features_unique)
  final_features_freq <- rep.int(0, length(final_features_unique))
  names(final_features_freq) <- final_features_unique
  for (population in Pops$final_features) {
    for (feature in population) {
      f.ind <- which(final_features_unique == feature) 
      final_features_freq[f.ind] <- final_features_freq[f.ind] + 1
    }
  }
  final_features.names <- vector()
  count <- 0
  while (length(final_features.names) < 40) {
    final_features.names <- c(final_features.names, final_features_unique[which(final_features_freq == length(Pops$final_features)-count)])
	count<-count+1
  }
  final_features.length <- length(final_features.names)
  #
  features_in_G1_unique <- unique(features_in_G1_unique)
  features_in_G1_freq <- rep.int(0, length(features_in_G1_unique))
  names(features_in_G1_freq) <- features_in_G1_unique
  for (population in Pops$features_in_G1) {
    for (feature in population) {
      f.ind <- which(features_in_G1_unique == feature)
      features_in_G1_freq[f.ind] <- features_in_G1_freq[f.ind] + 1
    }
  }
  best_features_in_G1.names <- vector()
  count <- 0
  while (length(best_features_in_G1.names) < 40) {
    best_features_in_G1.names <- c(best_features_in_G1.names, features_in_G1_unique[which(features_in_G1_freq == length(Pops$features_in_G1)-count)])
	count<-count+1
  }
  best_features_in_G1.length <- length(best_features_in_G1.names)
  #
  final_fitness.best <- mean(Pops$final_fitness)
  
  ### Final features rate
  correct_final_features.length <- length(which(final_features.names %in% DATA$major_snp_1))
  correct_final_features.rate <- correct_final_features.length / length(final_features.names)
  cat('\n correct_final_features.rate =', correct_final_features.rate, '\n')
  error_final_features.length <- length(final_features.names) - correct_final_features.length
  error_final_features.rate <- error_final_features.length / length(final_features.names)
  cat('error_final_features.rate =', error_final_features.rate, '\n')
  caught_final_features.length <- correct_final_features.length
  caught_final_features.rate <- caught_final_features.length / length(DATA$major_snp_1)
  cat('caught_final_features.rate =', caught_final_features.rate, '\n')
  lost_final_features.length <- length(DATA$major_snp_1) - caught_final_features.length
  lost_final_features.rate <- lost_final_features.length / length(DATA$major_snp_1)
  cat('lost_final_features.rate =', lost_final_features.rate, '\n')
  final_features.MR <- matrix(data = c(correct_final_features.length, lost_final_features.length, error_final_features.length,
                                       (ncol(DATA$m.probe) - correct_final_features.length - lost_final_features.length - error_final_features.length)), ncol = 2)
  colnames(final_features.MR) <- c("Major", "Rest")
  rownames(final_features.MR) <- c("Major", "Rest")
  final_heat.MR <- heatmaply(final_features.MR, dendrogram = "none", hide_colorbar = T, draw_cellnote = T, colors = rainbow(100, alpha = 0.75, start = 0, end = 0.15, rev = T),
                             cellnote_size = 22, grid_color = "white", column_text_angle = 0, xlab = 'Actual SNP', ylab = 'Predicted SNP',
                             main = 'Final features    .', fontsize_row = 12, fontsize_col = 12, subplot_widths = 0.4, subplot_heights = 1, 
                             cellnote_textposition = "middle center", cellnote_color = "black", colorbar_xpos = 3, 
                             file = paste0("results/", AnalyseName, "_final_heat.MR.html"))
  
  ### Best features in G1 rate
  correct_features_in_G1.length <- length(which(best_features_in_G1.names %in% DATA$major_snp_1))
  correct_features_in_G1.rate <- correct_features_in_G1.length / length(best_features_in_G1.names)
  cat('\n correct_features_in_G1.rate =', correct_features_in_G1.rate, '\n')
  error_features_in_G1.length <- length(best_features_in_G1.names) - correct_features_in_G1.length
  error_features_in_G1.rate <- error_features_in_G1.length / length(best_features_in_G1.names)
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
  BDE_Result <- as.data.frame(c(final_fitness.best, final_features.length, correct_final_features.rate, error_final_features.rate, caught_final_features.rate, 
                                lost_final_features.rate, best_features_in_G1.length, correct_features_in_G1.rate, error_features_in_G1.rate, 
                                caught_features_in_G1.rate, lost_features_in_G1.rate, BDE_time_print, 'BDE Parameters:', BDE_Parameters, 
                                'OBJ_FUNC Parameters:', OBJFUNC_Parameters), optional = T)
  row.names(BDE_Result) <- c('final_fitness.best', 'final_features.length', 'correct_final_features.rate', 
                             'error_final_features.rate','caught_final_features.rate', 'lost_final_features.rate', 'best_features_in_G1.length', 
                             'correct_features_in_G1.rate', 'error_features_in_G1.rate', 'caught_features_in_G1.rate', 
                             'lost_features_in_G1.rate', 'BDE_time', '###', names(BDE_Parameters), '####', names(OBJFUNC_Parameters))
  write.table(BDE_Result, file = paste0("results/", AnalyseName, "_BDE_Result.txt"), col.names = F)
  return(list(final_heat.MR=final_heat.MR, best_in_G1_heat.MR=best_in_G1_heat.MR))
}
