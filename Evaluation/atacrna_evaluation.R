library(Seurat)
## load data
dta_harmony_together<- readRDS("R_Code/results/harmony_results/atacrna/harmony_together.rds")
dta_harmony_seperate <- readRDS("R_Code/results/harmony_results/atacrna/harmony_seperate.rds")
dta_seurat <- readRDS("R_Code/results/seurat_results/atacrna/seurat_atacrna.rds")

batch <- c(rep("atac", 1047), rep("rna", 1047))
library(irlba)
library(Matrix)
dta_seurat_svd <- irlba(dta_seurat, nv = 100)
dta_seurat <- dta_seurat %*% dta_seurat_svd$v

## evaluation
source("Evaluation/utility.R")
# kBET
dta_harmony_toge_kBET <- kBET_fun(data = dta_harmony_together, batch = batch)
dta_harmony_sepe_kBET <- kBET_fun(data = dta_harmony_seperate, batch = batch)
dta_seurat_kBET <- kBET_fun(data = dta_seurat, batch = batch)
dta <- data.frame("reject rate" = c(dta_harmony_toge_kBET, dta_harmony_sepe_kBET, dta_seurat_kBET), 
                     "percentage of sample size" = rep(seq(0.05, 0.25, by = 0.05), 3),
                     "data type" = c(rep("harmony-together", 5), rep("harmony-seperate", 5), rep("seurat", 5)))
ggplot(data = dta, aes(x = percentage.of.sample.size, y = reject.rate)) +
  geom_line(aes(color = data.type)) +
  geom_point(aes(color = data.type))

# lisi
library(lisi)
batch <- as.data.frame(batch)
colnames(batch) <- "batch"
mean(compute_lisi(dta_harmony_together, meta_data = batch, label_colnames = "batch")[, 1])
mean(compute_lisi(dta_harmony_seperate, meta_data = batch, label_colnames = "batch")[, 1])
mean(compute_lisi(dta_seurat, meta_data = batch, label_colnames = "batch")[, 1])

# ARI
library(mclust)
library(stats)

dta_harmony_toge_kmeans <- kmeans(dta_harmony_together, nstart = 10, centers = 2)
adjustedRandIndex(dta_harmony_toge_kmeans$cluster, batch$batch)
dta_harmony_sepe_kmeans <- kmeans(dta_harmony_seperate, nstart = 10, centers = 2)
adjustedRandIndex(dta_harmony_sepe_kmeans$cluster, batch$batch)
dta_seurat_kmeans <- kmeans(dta_seurat, nstart = 10, centers = 2)
adjustedRandIndex(dta_seurat_kmeans$cluster, batch$batch)



