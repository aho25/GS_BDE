### Init markers ###
library(SoyNAM)
library(rrBLUP)

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

### Initialize data
n <- nrow(Markers)
p <- ncol(Markers)
probedata <- sample(1:n, 0.8 * n)
validata <- setdiff(1:n, probedata)
p.valid <- as.matrix(Pheno[validata,]$yield)
m.valid <- Markers[validata,]
p.probe <- as.matrix(Pheno[probedata,]$yield)
m.probe <- Markers[probedata,]
# 
### Check
set.seed(12)
m.rows = sample(n,300)
m.cols = sample(p,300)
m.probe = Markers[m.rows,m.cols]
p.probe = as.matrix(Pheno[m.rows,]$yield)
MARKERS = as.matrix(m.probe)

### make markers == -1, 0, 1
OFFSET <- 0.1
for (j in 1:ncol(MARKERS)) {
  for (i in 1:nrow(MARKERS)) {
    if (MARKERS[i,j] < -OFFSET) {
      MARKERS[i,j] <- -1
    } else if (MARKERS[i,j] > OFFSET) {
      MARKERS[i,j] <- 1
    } else {
      MARKERS[i,j] <- 0
    }
  }
}

### PHEN0 ###
# normal_distr <- rnorm(10^5)
# quantile_0.5 <- sapply(1:10^5, function(i) {
#   if (normal_distr[i] < 0.5) return(1)
#   else return(0)
# })
# pnorma_viborochnoe = sum(quantile_0.5)/10^5
# 
# qnorm(0.1)
# 

x1 <- runif(300, -3, 3)
x1 <- sort(x1)

weight_1 <- dnorm(x1)
weight_1.1 <- dnorm(x1, 0, 0.5)

priznak_1 <- MARKERS %*% weight_1
#priznak_1.1 <- MARKERS %*% weight_1.1

pheno_1 <- sapply(priznak_1, function(i) {
  (max(priznak_1) - min(priznak_1)) / (max(priznak_1) + 1 - i)
})

#pheno_1.1 <- sapply(priznak_1.1, function(i) {
 # (max(priznak_1.1) - min(priznak_1.1)) / (max(priznak_1.1) + 1 - i)
#})

#sum(sort(priznak_2)[150:151])/2

### multi trait ###
x2 <- runif(300, -3, 1)
x2 <- sort(x2)
weight_2 <- dnorm(x2, 1, 1.5)
priznak_2 <- MARKERS %*% weight_2
pheno_2 <- sapply(priznak_2, function(i) {
  (max(priznak_2) - min(priznak_2)) / (max(priznak_2) + 1 - i)
})

cor(pheno_1, pheno_2)
# Make pheno_2 correlated to pheno_1

ro <- 1 #cor_coef
x2.1 <- ro*x1 + sqrt(1-ro^2)*x2

weight_2.1 <- dnorm(x2.1, 1, 1.5)
priznak_2.1 <- MARKERS %*% weight_2.1
pheno_2.1 <- sapply(priznak_2.1, function(i) {
  (max(priznak_2.1) - min(priznak_2.1)) / (max(priznak_2.1) + 1 - i)
})


### Concanate pheno
PHENO <- cbind(pheno_1, pheno_2.1)

n <- 1:nrow(MARKERS)
probedata <- sample(1:n, 0.8)
validata <- setdiff(1:n, probedata)
m.train <- MARKERS[probedata,]
p.train <- PHENO[probedata,]
m.test <- MARKERS[validata,]
p.test <- PHENO[validata,]

