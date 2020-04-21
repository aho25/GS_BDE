library(SoyNAM)

Soy_init_MT <- function() {
  
  ### Load data
  data(soybase)
  dim(gen.qa)
  dim(data.line.qa)
  
  ### Data imputation
  idx_arr <- rownames(gen.qa) %in% data.line.qa$strain
  Markers.src = gen.qa[idx_arr,]
  Markers.src[Markers.src == 0] <- -1
  Markers.src[Markers.src == 1] <- 0
  Markers.src[Markers.src == 2] <- 1
  summary(Markers.src[,1:10])
  sum(is.na(Markers.src))/nrow(Markers.src)/ncol(Markers.src)
  impute <- A.mat(Markers.src, min.MAF = 0.05, max.missing = 0.5, impute.method = "mean", return.imputed = T)
  markers_imputed <- as.data.frame(impute$imputed)
  dim(markers_imputed)
  summary(markers_imputed[,1:10])
  
  ### Order data
  idx_arr <- data.line.qa$strain %in% rownames(markers_imputed)
  Pheno.src <- data.line.qa[idx_arr,]
  dim(Pheno.src)
  podmnogestvo <- unique(c(which(is.na(Pheno.src$yield)), which(is.na(Pheno.src$moisture)), which(is.na(Pheno.src$protein)),
                           which(is.na(Pheno.src$oil)), which(is.na(Pheno.src$fiber)), which(is.na(Pheno.src$size))))  #delete 'NA' from traits
  Pheno.src <- Pheno.src[-podmnogestvo,]
  dim(Pheno.src)
  idx_vec <- match(Pheno.src$strain, rownames(markers_imputed))
  Markers.src <- markers_imputed[idx_vec,]
  dim(Markers.src)
  
  ### Location: 'NE', Year: 2012
  USED_YEAR = c(2012)
  podmnogestvo <- Pheno.src$year %in% USED_YEAR
  Pheno <- Pheno.src[podmnogestvo,]
  Markers <- Markers.src[podmnogestvo,]
  USED_LOCATION = c('NE')
  podmnogestvo <- Pheno$location %in% USED_LOCATION
  Pheno <- Pheno[podmnogestvo,]
  Markers <- Markers[podmnogestvo,]
  dim(Markers)
  dim(Pheno)
  
  Markers.2012.NE <- Markers
  Pheno.2012.NE <- Pheno
  
  ### Check correlation rates between traits
  Pheno_traits <- cbind(Pheno$yield, Pheno$moisture, Pheno$protein, Pheno$oil, Pheno$fiber, Pheno$size)
  colnames(Pheno_traits) <- c('yield', 'moisture', 'protein', 'oil', 'fiber', 'size')
  rownames(Pheno_traits) <- Pheno$strain
  cor(Pheno_traits)
  Pheno_corTraits <- Pheno_traits[,-c(2)] #leave traits correlated with yield
  
  ### Initialize data
  n <- nrow(Markers)
  p <- ncol(Markers)
  probedata <- sample(1:n, 0.8 * n)
  validata <- setdiff(1:n, probedata)
  p.valid <- Pheno_corTraits[validata,]
  m.valid <- Markers[validata,]
  p.probe <- Pheno_corTraits[probedata,]
  m.probe <- Markers[probedata,]
  
  # ### Check
  # set.seed(12)
  # m.rows = sample(n,300)
  # m.cols = sample(p,300)
  # m.probe = Markers[m.rows,m.cols]
  # p.probe = Pheno_corTraits[m.rows,]
  # MARKERS = m.probe
  # PHENO = p.probe
  
  Init_data <- list(p.valid=p.valid, m.valid=m.valid, p.probe=p.probe, m.probe=m.probe)
  return(Init_data)
}
