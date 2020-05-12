Synthetic_data_ST2 <- function() {
  
  ### Load MARKERS ###
  MARKERS <- read.csv('data_csv/markers_synthetic.csv')
  rownames(MARKERS) <- MARKERS[,1]
  MARKERS <- MARKERS[,-1]
  
  ### Load PHENO ###
  PHENO <- read.csv('data_csv/pheno2_synthetic_st_001.csv')#data to change
  PHENO <- PHENO[,-1]
  names(PHENO) <- rownames(MARKERS)
  
  major_snp_1 <- read.csv('data_csv/major_snp1_synthetic_st_001.csv')#change here too
  major_snp_1 <- as.vector(major_snp_1[,-1])
  
  MARKERS <- as.data.frame(MARKERS)
  PHENO <- as.matrix(PHENO)

  Init_data <- list(major_snp_1=major_snp_1, p.probe=PHENO, m.probe=MARKERS)
  return(Init_data)
}