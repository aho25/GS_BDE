#data_ST2()
Synthetic_data_ST <- function(k) {

  ### Load MARKERS ###
  MARKERS <- read.csv('data_csv/synthetic/markers_synthetic.csv')
  rownames(MARKERS) <- MARKERS[,1]
  MARKERS <- MARKERS[,-1]
  
  ### Load PHENO ###
  PHENO <- read.csv(paste0('data_csv/synthetic/pheno_synthetic_st_', k, '.csv'))
  PHENO <- PHENO[,-1]
  names(PHENO) <- rownames(MARKERS)
  
  weight_1 <- read.csv(paste0('data_csv/synthetic/weight_1_synthetic_st_', k, '.csv'))
  weight_1 <- weight_1[,-1]
  names(weight_1) <- colnames(MARKERS)
  
  major_snp_1 <- read.csv(paste0('data_csv/synthetic/major_snp_1_synthetic_st_', k, '.csv'))
  major_snp_1 <- as.vector(major_snp_1[,-1])
  
  MARKERS <- as.data.frame(MARKERS)
  PHENO <- as.matrix(PHENO)

  Init_data <- list(p.probe=PHENO, m.probe=MARKERS, weight_1=weight_1, major_snp_1=major_snp_1)
  return(Init_data)
}