library(SoyNAM)
library(rrBLUP)
#можно изменять ср и дисперсии
Synthetic_data_MT <- function(OFFSET) {
  
  ### Init markers ###
  data(soybase)
  
<<<<<<< HEAD
  ### Load PHENO ###
  PHENO <- read.csv('data_csv/pheno_synthetic_mt_001.csv')#data to change
  PHENO <- PHENO[,-1]
  row.names(PHENO) <- row.names(MARKERS)
  
  major_snp_1 <- read.csv('data_csv/major_snp1_synthetic_mt_001.csv')#data to change
  major_snp_1 <- as.vector(major_snp_1[,-1])
=======
  set.seed(12)
  Markers <- gen.qa
  n <- nrow(Markers)
  p <- ncol(Markers)
  m.rows = sample(n, 500)
  m.cols = sample(p, 500)
  Markers = as.matrix(Markers[m.rows, m.cols])
  
  ### Data imputation
  Markers[Markers == 0] <- -1
  Markers[Markers == 1] <- 0
  Markers[Markers == 2] <- 1
  impute <-
    A.mat(
      Markers,
      min.MAF = 0.05,
      max.missing = 0.5,
      impute.method = "mean",
      return.imputed = T
    )
  MARKERS <- as.data.frame(impute$imputed)
  dim(MARKERS)
  
  ### Truncate data
  n <- nrow(MARKERS)
  p <- ncol(MARKERS)
  m.rows = sample(n, 300)
  m.cols = sample(p, 300)
  MARKERS = as.matrix(MARKERS[m.rows, m.cols])
  
  ### Make MARKERS == -1, 0, 1
  for (j in 1:ncol(MARKERS)) {
    for (i in 1:nrow(MARKERS)) {
      if (MARKERS[i, j] < -OFFSET) {
        MARKERS[i, j] <- -1
      } else if (MARKERS[i, j] > OFFSET) {
        MARKERS[i, j] <- 1
      } else {
        MARKERS[i, j] <- 0
      }
    }
  }
  
  ### Create PHENO ###
  set.seed(12)
  major_len_1 <- 20
  major_ind_1 <- sample(length(m.cols), major_len_1)
  weight_1 <- rnorm(length(m.cols), 3, 2)
  weight_1[major_ind_1] <- weight_1[major_ind_1] + rnorm(major_len_1, 15, 2)
  names(weight_1) <- colnames(MARKERS)
  major_snp_1 <- names(sort(weight_1[major_ind_1]))
  priznak_1 <- MARKERS %*% weight_1
  pheno_1 <- priznak_1 + 2*abs(min(priznak_1)) + rnorm(m.rows, 0, 10)
  # max(pheno_1)
  # min(pheno_1)
  # mean(pheno_1)
  # scale_1 <- 1
  # pheno_1 <- pheno_1 * scale_1
  
  ### Multi trait
  major_len_2 <- 40
  major_ind_2 <- sample(length(m.cols), major_len_2)
  weight_2 <- rnorm(length(m.cols), 2, 2)
  weight_2[major_ind_2] <- weight_2[major_ind_2] + rnorm(major_len_2, 25, 4)
  names(weight_2) <- colnames(MARKERS)
  major_snp_2 <- names(sort(weight_2[major_ind_2]))
  priznak_2 <- MARKERS %*% weight_2
  pheno_2 <- priznak_2 + 2*abs(min(priznak_2)) + rnorm(m.rows, 0, 10)
  # max(pheno_2)
  # min(pheno_2)
  # mean(pheno_2)
  scale_2 <- 0.05
  pheno_2 <- pheno_2 * scale_2
  #cor(pheno_1, pheno_2) #выяснить при каком значении корелляции, улучшается/ухудшается работа алгоритма МТ. 
  
  ### Pheno_3
  major_len_3 <- 10
  major_ind_3 <- sample(length(m.cols), major_len_3)
  weight_3 <- rnorm(length(m.cols), 1, 3)
  weight_3[major_ind_3] <- weight_3[major_ind_3] + rnorm(major_len_3, 30, 5)
  names(weight_3) <- colnames(MARKERS)
  major_snp_3 <- names(sort(weight_3[major_ind_3]))
  priznak_3 <- MARKERS %*% weight_3
  pheno_3 <- priznak_3 + 2*abs(min(priznak_3)) + rnorm(m.rows, 0, 10)
  #max(pheno_3)
  #min(pheno_3)
  #mean(pheno_3)
  scale_3 <- 0.2
  pheno_3 <- pheno_3 * scale_3
  #cor(pheno_1, pheno_3)
  
  ### Concanate pheno
  PHENO <- cbind(pheno_1, pheno_2, pheno_3)
  colnames(PHENO) <- c('yield', 'pheno_2', 'pheno_3')
>>>>>>> parent of 33a3ed4... Syn data update, complete data sets + results
  
  MARKERS <- as.data.frame(MARKERS)
  PHENO <- as.matrix(PHENO)
  Init_data <- list(major_snp_1=major_snp_1, major_snp_2=major_snp_2, major_snp_3=major_snp_3, p.probe=PHENO, m.probe=MARKERS)
  return(Init_data)
}