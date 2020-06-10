# data: a matrix (rows: samples, columns: features (genes))
# k0: number of nearest neighborhood
kBET_fun <- function(data, batch, k0 = NULL, subset_size = 1, seed.use = 1234){
  library(kBET)
  library(FNN)
  if (is.null(k0) == TRUE){
    k0_max <- nrow(data)
    k0 <- seq(0.05, 0.25, by = 0.05) * k0_max
    k0 <- floor(k0)
  }
  set.seed(seed.use)
  subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
  data <- data[subset_id,,drop = FALSE]
  batch <- batch[subset_id]
  
  reject_rate <- rep(NA, length(k0))
  for (i in 1:length(k0)) {
    k <- k0[i]
    knn <- get.knn(data, k=k, algorithm = 'cover_tree')
    batch.estimate <- kBET(data, batch, k = k, knn = knn, plot = FALSE)
    reject_rate[i] <- mean(batch.estimate$stats$kBET.observed)
  }
  return(reject_rate)
}





