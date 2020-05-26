#pheno2_syntetic
### Set noise coefficient
k <- 0.1

### Download markers
MARKERS <- read.csv('data_csv/synthetic/markers_synthetic.csv')
rownames(MARKERS) <- MARKERS[,1]
MARKERS <- MARKERS[,-1]

### Create PHENO ###
set.seed(12)
major_len_1 <- 20
major_ind_1 <- sample(ncol(MARKERS), major_len_1)
weight_1 <- rep(0,ncol(MARKERS))
weight_1[major_ind_1] <- weight_1[major_ind_1] + rnorm(major_len_1, 10, 1)
names(weight_1) <- colnames(MARKERS)
major_snp_1 <- names(sort(weight_1[major_ind_1], decreasing = T))
priznak_1 <- as.vector(crossprod(t(MARKERS),weight_1))
PHENO <- priznak_1
PHENO <- PHENO + rnorm(nrow(MARKERS), 0, 1)*k

# max(PHENO)
# min(PHENO)
# mean(PHENO)
# scale_1 <- 1
# PHENO <- PHENO * scale_1

write.csv(PHENO, file = paste0('data_csv/synthetic/pheno_synthetic_st_', k, '.csv'))
write.csv(weight_1, file = paste0('data_csv/synthetic/weight_1_synthetic_st_', k, '.csv'))
write.csv(major_snp_1, file = paste0('data_csv/synthetic/major_snp_1_synthetic_st_', k, '.csv'))

### Control_discrete
# #predict marker effects
# ans <- mixed.solve(PHENO,Z=MARKERS)
# best_features <- rownames(sort(ans$u, decreasing = T))
# best_features_len <- length(which(best_features[1:20] %in% DATA$major_snp_1))
