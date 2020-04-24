library(SoyNAM)
library(rrBLUP)

Synthetic_data_MT <- function(OFFSET) {
  
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
  x1 <- sort(runif(m.cols,-3, 3))
  weight_1 <- dnorm(x1)
  priznak_1 <- MARKERS %*% weight_1
  pheno_1 <- sapply(priznak_1, function(i) {
    (max(priznak_1) - min(priznak_1)) / (max(priznak_1) + 1 - i)
  })
  
  ### Multi trait
  x2 <- sort(runif(m.cols,-3, 1))
  weight_2 <- dnorm(x2, 1, 1.5)
  priznak_2 <- MARKERS %*% weight_2
  pheno_2 <- sapply(priznak_2, function(i) {
    (max(priznak_2) - min(priznak_2)) / (max(priznak_2) + 1 - i)
  })
  #cor(pheno_1, pheno_2)
  
  ### Make pheno_2 correlated to pheno_1
  ro_2 <- 0.3 #cor_coef
  pheno_2 <- ro_2 * pheno_1 + sqrt(1 - ro_2 ^ 2) * pheno_2
  #cor(pheno_1, pheno_2)
  
  ### Pheno_3
  x3 <- sort(runif(m.cols, -2, 4))
  weight_3 <- dnorm(x3, -1, 1.5)
  priznak_3 <- MARKERS %*% weight_3
  pheno_3 <- sapply(priznak_2, function(i) {
    (max(priznak_2) - min(priznak_2)) / (max(priznak_2) + 1 - i)
  })
  #cor(pheno_1, pheno_3)
  ro_3 <- 0.2 #cor_coef
  pheno_3 <- ro_3 * pheno_1 + sqrt(1 - ro_3 ^ 2) * pheno_3
  #cor(pheno_1, pheno_3)
  
  ### Concanate pheno
  PHENO <- cbind(pheno_1, pheno_2, pheno_3)
  colnames(PHENO) <- c('yield', 'pheno_2', 'pheno_3')
  
  n <- nrow(MARKERS)
  MARKERS <- as.data.frame(MARKERS)
  PHENO <- as.matrix(PHENO)
  probedata <- sample(1:n, 0.8 * n)
  validata <- setdiff(1:n, probedata)
  m.probe <- MARKERS[probedata, ]
  p.probe <- as.matrix(PHENO[probedata, ])
  m.valid <- MARKERS[validata, ]
  p.valid <- as.matrix(PHENO[validata, ])
  Init_data <- list(p.valid=p.valid, m.valid=m.valid, p.probe=p.probe, m.probe=m.probe)
  return(Init_data)
}