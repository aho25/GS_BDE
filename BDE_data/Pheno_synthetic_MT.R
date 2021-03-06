### Set noise coefficient
k <- 1

### Download markers
MARKERS <- read.csv('data_csv/synthetic/markers_synthetic.csv')
rownames(MARKERS) <- MARKERS[,1]
MARKERS <- MARKERS[,-1]

### Create PHENO ###
set.seed(12)
major_len_1 <- 20
major_ind_1 <- sample(ncol(MARKERS), major_len_1)
weight_1 <- rnorm(ncol(MARKERS), 0, 0.001)
weight_1[major_ind_1] <- weight_1[major_ind_1] + rnorm(major_len_1, 10, 1)
names(weight_1) <- colnames(MARKERS)
major_snp_1 <- names(sort(weight_1[major_ind_1], decreasing = T))
priznak_1 <- as.vector(crossprod(t(MARKERS),weight_1))
PHENO1 <- priznak_1
PHENO1 <- PHENO1 + rnorm(nrow(MARKERS), 0, 1)*sd(PHENO1)*k
# max(PHENO1)
# min(PHENO1)
# mean(PHENO1)
# scale_1 <- 1
# PHENO1 <- PHENO1 * scale_1

write.csv(weight_1, file = paste0('data_csv/synthetic/weight_1_synthetic_mt_', k, '.csv'))
write.csv(major_snp_1, file = paste0('data_csv/synthetic/major_snp_1_synthetic_mt_', k, '.csv'))

### Create PHENO2 ###
major_len_2 <- 20
major_len_2_in_major_len_1 <- 15
major_ind_2 <- sample(ncol(MARKERS), major_len_2)
major_ind_2 <- major_ind_2[-which(major_ind_2 %in% major_ind_1)]
major_ind_2 <- c(major_ind_1[sample(major_len_1, major_len_2_in_major_len_1)], 
                 major_ind_2[sample(length(major_ind_2), major_len_2 - major_len_2_in_major_len_1)])
weight_2 <- rnorm(ncol(MARKERS), 0, 0.001)
weight_2[major_ind_2] <- weight_2[major_ind_2] + rnorm(major_len_2, 15, 3)
names(weight_2) <- colnames(MARKERS)
major_snp_2 <- names(sort(weight_2[major_ind_2]))
priznak_2 <- as.vector(crossprod(t(MARKERS),weight_2))
PHENO2 <- priznak_2
PHENO2 <- PHENO2 + rnorm(nrow(MARKERS), 0, 1)*sd(PHENO2)*k
# max(PHENO2)
# min(PHENO2)
# mean(PHENO2)
scale_2 <- 0.5
PHENO2 <- PHENO2 * scale_2

write.csv(weight_2, file = paste0('data_csv/synthetic/weight_2_synthetic_mt_', k, '.csv'))

### Create PHENO3 ###
major_len_3 <- 25
major_len_3_in_major_len_1 <- 10
major_ind_3 <- sample(ncol(MARKERS), major_len_3)
major_ind_3 <- major_ind_3[-which(major_ind_3 %in% major_ind_1)]
major_ind_3 <- c(major_ind_1[sample(major_len_1, major_len_3_in_major_len_1)], 
                 major_ind_3[sample(length(major_ind_3), major_len_3 - major_len_3_in_major_len_1)])
weight_3 <- rnorm(ncol(MARKERS), 0, 0.001)
weight_3[major_ind_3] <- weight_3[major_ind_3] + rnorm(major_len_3, 10, 2)
names(weight_3) <- colnames(MARKERS)
major_snp_3 <- names(sort(weight_3[major_ind_3]))
priznak_3 <- as.vector(crossprod(t(MARKERS),weight_3))
PHENO3 <- priznak_3
PHENO3 <- PHENO3 + rnorm(nrow(MARKERS), 0, 1)*sd(PHENO2)*k
# max(PHENO3)
# min(PHENO3)
# mean(PHENO3)
scale_3 <- 0.1
PHENO3 <- PHENO3 * scale_3

# ro_3 <- 0.7
# weight_3 <- ro_3 * weight_1 + (1 - ro_3) * weight_3
# cor(weight_1, weight_3)

write.csv(weight_3, file = paste0('data_csv/synthetic/weight_3_synthetic_mt_', k, '.csv'))

### Concatanate pheno ###
PHENO <- cbind(PHENO1, PHENO2, PHENO3)
colnames(PHENO) <- c('yield', 'pheno_2', 'pheno_3')
row.names(PHENO) <- row.names(MARKERS)
write.csv(PHENO, file = paste0('data_csv/synthetic/pheno_synthetic_mt_', k, '.csv'))

### Create weight-cov-matrix ###
weight_cov <- cov(cbind(weight_1, weight_2, weight_3))
write.csv(weight_cov, file = paste0('data_csv/synthetic/weight_cov_mt_', k, '.csv'))
