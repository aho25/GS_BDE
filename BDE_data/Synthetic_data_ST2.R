Synthetic_data_ST2 <- function(k) {
  
  ### Load MARKERS ###
  MARKERS <- read.csv('data_csv/synthetic/markers_synthetic.csv')
  rownames(MARKERS) <- MARKERS[,1]
  MARKERS <- MARKERS[,-1]
  
  ### Load PHENO ###
  PHENO <- read.csv(paste0('data_csv/synthetic/pheno_synthetic_st_', k, '.csv'))
  PHENO <- PHENO[,-1]
  names(PHENO) <- rownames(MARKERS)
  index<-sample(nrow(MARKERS),50)
  MARKERS_VALID<-MARKERS[index,]  
  PHENO_VALID<-PHENO[index]  
  MARKERS_TRAIN<-MARKERS[-index, ]
  PHENO_TRAIN<-PHENO[-index]
  MARKERS_VALID <- as.data.frame(MARKERS)
  PHENO_VALID <- as.matrix(PHENO)
  MARKERS_TRAIN <- as.data.frame(MARKERS)
  PHENO_TRAIN <- as.matrix(PHENO)
  Init_data <- list(p.valid=PHENO_VALID, m.valid=MARKERS_VALID, 
                    p.probe=PHENO_TRAIN, m.probe=MARKERS_TRAIN)
  return(Init_data)
}
