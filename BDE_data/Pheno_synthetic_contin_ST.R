##pheno_synthetic генерит данные и веса из нормального распределения
### Download markers
MARKERS <- read.csv('data_csv/synthetic/markers_synthetic.csv')
rownames(MARKERS) <- MARKERS[,1]
MARKERS <- MARKERS[,-1]

### Create PHENO ###
set.seed(12)
#random phenotypes
weight_1 <- rnorm(400)
g <- as.vector(crossprod(t(MARKERS),weight_1))
h2 <- 1  #heritability
PHENO <- g + rnorm(400,mean=0,sd=sqrt((1-h2)/h2*var(g)))
names(PHENO) <- rownames(MARKERS)
names(weight_1) <- colnames(MARKERS)

write.csv(PHENO, file = 'data_csv/synthetic/pheno_synthetic_contin_st.csv')
write.csv(weight_1, file = 'data_csv/synthetic/weight_1_synthetic_contin_st.csv')
