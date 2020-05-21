library(SoyNAM)
library(rrBLUP)

### Init markers ###
data(soybase)
OFFSET <- 0.1

set.seed(12)
Markers <- gen.qa
n <- nrow(Markers)
p <- ncol(Markers)
m.rows = sample(n, 500)
m.cols = sample(p, 800)
Markers = as.matrix(Markers[m.rows, m.cols])

### Data imputation
Markers[Markers == 0] <- -1
Markers[Markers == 1] <- 0
Markers[Markers == 2] <- 1
impute <-
  A.mat(
    Markers,
    min.MAF = 0.05,
    max.missing = 0.5,
    impute.method = "mean",
    return.imputed = T
  )
MARKERS <- as.data.frame(impute$imputed)
dim(MARKERS)

### Truncate data
n <- nrow(MARKERS)
p <- ncol(MARKERS)
m.rows = sample(n, 400)
m.cols = sample(p, 400)
MARKERS = as.matrix(MARKERS[m.rows, m.cols])

### Make MARKERS == -1, 0, 1
for (j in 1:ncol(MARKERS)) {
  for (i in 1:nrow(MARKERS)) {
    if (MARKERS[i, j] < -OFFSET) {
      MARKERS[i, j] <- -1
    } else if (MARKERS[i, j] > OFFSET) {
      MARKERS[i, j] <- 1
    } else {
      MARKERS[i, j] <- 0
    }
  }
}
write.csv(MARKERS, file = 'data_csv/synthetic/markers_synthetic.csv')
