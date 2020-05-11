library(SoyNAM)
library(rrBLUP)

Synthetic_data_ST <- function(OFFSET) {
  
  ### Init markers ###
  data(soybase)
  
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
  major_len_1 <- 30  #количество главных фич(влияют на результат(вариация по синтетическим данным))
  major_ind_1 <- sample(length(m.cols), major_len_1)
  weight_1 <- rnorm(length(m.cols), 0, 1)#from normal distribution, select random value 
  weight_1[major_ind_1] <- weight_1[major_ind_1] #add normal distr to weight
  names(weight_1) <- colnames(MARKERS)
  major_snp_1 <- names(sort(weight_1[major_ind_1]))
  priznak_1 <- MARKERS %*% weight_1
  PHENO <- priznak_1 + 2*abs(min(priznak_1))
  PHENO <- PHENO + runif(m.rows)*rnorm(m.rows,mean(PHENO)/10, 1)
  max(PHENO)
  min(PHENO)
  mean(PHENO)
  scale_1 <- 1
  PHENO <- PHENO * scale_1
  
  MARKERS <- as.data.frame(MARKERS)
  PHENO <- as.matrix(PHENO)
  Init_data <- list(major_snp_1=major_snp_1, p.probe=PHENO, m.probe=MARKERS)
  return(Init_data)
}
