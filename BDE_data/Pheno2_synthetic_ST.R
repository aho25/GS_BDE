### Download markers
MARKERS <- read.csv('data_csv/markers_synthetic.csv')
rownames(MARKERS) <- MARKERS[,1]
MARKERS <- MARKERS[,-1]

### Create PHENO ###
set.seed(12)
major_len_1 <- 20
major_ind_1 <- sample(ncol(MARKERS), major_len_1)
weight_1 <- rnorm(ncol(MARKERS), 0, 0.1)
weight_1[major_ind_1] <- weight_1[major_ind_1] + rnorm(major_len_1, 10, 1)#is variance can makes a difference to results?
names(weight_1) <- colnames(MARKERS)
major_snp_1 <- names(sort(weight_1[major_ind_1], decreasing = T))
priznak_1 <- as.vector(crossprod(t(MARKERS),weight_1))
PHENO <- priznak_1 + 2*abs(min(priznak_1))
PHENO <- PHENO + rnorm(nrow(MARKERS), 0, mean(PHENO)*0.1)#create datasets with different variance (0.01; 0.1 )
# max(PHENO)
# min(PHENO)
# mean(PHENO)
# scale_1 <- 1
# PHENO <- PHENO * scale_1

write.csv(PHENO, file = 'data_csv/pheno2_synthetic_st_01.csv')
write.csv(major_snp_1, file = 'data_csv/major_snp1_synthetic_st_01.csv')

# #predict marker effects
# ans <- mixed.solve(PHENO,Z=MARKERS)
# best_features <- rownames(sort(ans$u, decreasing = T))
# best_features_len <- length(which(best_features[1:20] %in% DATA$major_snp_1))

