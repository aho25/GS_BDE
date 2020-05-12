Soy_init_MT <- function() {
  
  ### Load data
  p.valid <- read.csv('data_csv/soy_mt_pheno_valid.csv')
  m.valid <- read.csv('data_csv/soy_mt_markers_valid.csv')
  p.valid <- as.matrix(p.valid[,-1])
  row.names(p.valid) <- m.valid[,1]
  m.valid <- m.valid[,-1]
  row.names(m.valid) <- row.names(p.valid)
  
  p.probe <- read.csv('data_csv/soy_mt_pheno_probe.csv')
  m.probe <- read.csv('data_csv/soy_mt_markers_probe.csv')
  p.probe <- as.matrix(p.probe[,-1])
  row.names(p.probe) <- m.probe[,1]
  m.probe <- m.probe[,-1]
  row.names(m.probe) <- row.names(p.probe)

  Init_data <- list(p.valid=p.valid, m.valid=m.valid, p.probe=p.probe, m.probe=m.probe)
  return(Init_data)
}
