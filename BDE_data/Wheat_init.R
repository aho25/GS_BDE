Wheat_init <- function() {
  
  ### Load data
  traits <- read.table('download_LRQW/traits.txt', header = T)
  snpfile <- read.table('download_LRQW/snpfile.txt')
  
  ### Data imputation
  impute <- A.mat(snpfile, min.MAF = 0.05, max.missing = 0.5, impute.method = "mean", return.imputed = T)
  markers_imputed <- as.data.frame(impute$imputed)
  
  ### Initialize data
  idx_arr <- traits$line %in% rownames(markers_imputed)
  pheno <- traits$thousand_kernel_weight[idx_arr]
  pheno <- pheno[!is.na(pheno)]
  names(pheno) <- traits$line[idx_arr]
  idx_vec <- match(names(pheno), rownames(markers_imputed))
  markers <- markers_imputed[idx_vec,]
  #summary(markers[,1:6])
  #any(is.na(markers))
  n <- nrow(markers)
  p <- ncol(markers)
  probedata <- sample(1:n, 0.8 * n)
  validata <- setdiff(1:n, probedata)
  p.valid <- as.matrix(pheno[validata])
  m.valid <- markers[validata,]
  p.probe <- as.matrix(pheno[probedata])
  m.probe <- markers[probedata,]
  
  # ### Check
  # set.seed(12)
  # m.rows = sample(n,300)
  # m.cols = sample(p,300)
  # p.probe = as.matrix(pheno[m.rows])
  # m.probe = markers[m.rows,m.cols]
  # PHENO = p.probe
  # MARKERS = m.probe
  
  Init_data <- list(p.valid=p.valid, m.valid=m.valid, p.probe=p.probe, m.probe=m.probe)
  return(Init_data)
}
