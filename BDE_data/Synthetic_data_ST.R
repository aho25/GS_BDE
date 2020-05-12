Synthetic_data_ST <- function() {
  
  ### Load MARKERS ###
  MARKERS <- read.csv('data_csv/markers_synthetic.csv')
  rownames(MARKERS) <- MARKERS[,1]
  MARKERS <- MARKERS[,-1]
  
  ### Load PHENO ###
  PHENO <- read.csv('data_csv/pheno_synthetic.csv') #rnorm data set
  PHENO <- PHENO[,-1]
  names(PHENO) <- rownames(MARKERS)
  
  weight_1 <- read.csv('data_csv/weight1_synthetic.csv')
  weight_1 <- weight_1[,-1]
  names(weight_1) <- colnames(MARKERS)
  
  MARKERS <- as.data.frame(MARKERS)
  PHENO <- as.matrix(PHENO)
  
  Init_data <- list(p.probe=PHENO, m.probe=MARKERS, weight=weight_1)
  return(Init_data)
}