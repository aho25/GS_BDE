Synthetic_data_MT <- function(k) {
  
  ### Load MARKERS ###
  MARKERS <- read.csv('data_csv/synthetic/markers_synthetic.csv')
  rownames(MARKERS) <- MARKERS[,1]
  MARKERS <- MARKERS[,-1]
  
  ### Load PHENO ###
  PHENO <- read.csv(paste0('data_csv/synthetic/pheno_synthetic_mt_', k, '.csv'))
  PHENO <- PHENO[,-1]
  row.names(PHENO) <- row.names(MARKERS)
  
  weight_1 <- read.csv(paste0('data_csv/synthetic/weight_1_synthetic_mt_', k, '.csv'))
  weight_1 <- weight_1[,-1]
  names(weight_1) <- colnames(MARKERS)
  
  major_snp_1 <- read.csv(paste0('data_csv/synthetic/major_snp_1_synthetic_mt_', k, '.csv'))
  major_snp_1 <- as.vector(major_snp_1[,-1])
  
  weight_cov <- read.csv(paste0('data_csv/synthetic/weight_cov_mt_', k, '.csv'))
  weight_cov <- as.vector(weight_cov[,-1])
  
  MARKERS <- as.data.frame(MARKERS)
  PHENO <- as.matrix(PHENO)
  Init_data <- list(p.probe=PHENO, m.probe=MARKERS, weight_1=weight_1, major_snp_1=major_snp_1, weight_cov=weight_cov)
  return(Init_data)
}