library(Seurat)
rm(list = ls())
## load data
dta_harmony<- readRDS("R_Code/results/harmony_results/pbmc/pbmc_harmony.rds")
metadata_harmony <- readRDS("R_Code/results/harmony_results/pbmc/meta_data.rds")
dta_seurat <- readRDS("R_Code/results/seurat_results/pbmc/pbmc_seurat.rds")
metadata_seurat <- readRDS("R_Code/results/seurat_results/pbmc/meta_data.rds")

## evaluation
source("Evaluation/utility.R")
# kBET
dta_harmony_kBET <- kBET_fun(data = dta_harmony, batch = metadata_harmony$batchlb, subset_size = 0.25)
dta_seurat_kBET <- kBET_fun(data = dta_seurat, batch = metadata_seurat$batchlb, subset_size = 0.25)
dta <- data.frame("reject rate" = c(dta_harmony_kBET, dta_seurat_kBET), 
                  "percentage of sample size" = rep(seq(0.05, 0.25, by = 0.05), 2),
                  "data type" = c(rep("harmony", 5), rep("seurat", 5)))
ggplot(data = dta, aes(x = percentage.of.sample.size, y = reject.rate)) +
  geom_line(aes(color = data.type)) +
  geom_point(aes(color = data.type))

# lisi
library(lisi)
# harmony
batch <- as.data.frame(metadata_harmony$batchlb)
colnames(batch) <- "batch"
celltype <- as.data.frame(metadata_harmony$cell_type)
colnames(celltype) <- "celltype"
batch_harmony <- mean(compute_lisi(dta_harmony, meta_data = batch, label_colnames = "batch")[, 1])
celltype_harmony <- mean(compute_lisi(dta_harmony, meta_data = celltype, label_colnames = "celltype")[, 1])
# seurat
batch <- as.data.frame(metadata_seurat$batchlb)
colnames(batch) <- "batch"
celltype <- as.data.frame(metadata_seurat$cell_type)
colnames(celltype) <- "celltype"
batch_seurat <- mean(compute_lisi(dta_seurat, meta_data = batch, label_colnames = "batch")[, 1])
celltype_seurat <- mean(compute_lisi(dta_seurat, meta_data = celltype, label_colnames = "celltype")[, 1])
# plot results
lisi_res <- data.frame("batch" = c(batch_harmony, batch_seurat), "celltype" = c(celltype_harmony, celltype_seurat), 
                       "datatype" = c("harmony", "seurat"))
ggplot(data = lisi_res, aes(x = batch, y = celltype)) +
  geom_point(aes(color = datatype))


# ARI
library(mclust)
library(stats)
# harmony
dta_harmony_kmeans_batch <- kmeans(dta_harmony, nstart = 10, centers = length(unique(metadata_harmony$batchlb)))
ari_batch_harmony <- adjustedRandIndex(dta_harmony_kmeans_batch$cluster, metadata_harmony$batchlb)
dta_harmony_kmeans_celltype <- kmeans(dta_harmony, nstart = 10, centers = length(unique(metadata_harmony$cell_type)))
ari_celltype_harmony <- adjustedRandIndex(dta_harmony_kmeans_celltype$cluster, metadata_harmony$cell_type)
# seurat
dta_seurat_kmeans_batch <- kmeans(dta_seurat, nstart = 10, centers = length(unique(metadata_seurat$batchlb)))
ari_batch_seurat <- adjustedRandIndex(dta_seurat_kmeans_batch$cluster, metadata_seurat$batchlb)
dta_seurat_kmeans_celltype <- kmeans(dta_seurat, nstart = 10, centers = length(unique(metadata_seurat$cell_type)))
ari_celltype_seurat <- adjustedRandIndex(dta_seurat_kmeans_celltype$cluster, metadata_seurat$cell_type)
# plot results
ari_res <- data.frame("batch" = c(ari_batch_harmony, ari_batch_seurat), "celltype" = c(ari_celltype_harmony, ari_celltype_seurat), 
                       "datatype" = c("harmony", "seurat"))
ggplot(data = ari_res, aes(x = batch, y = celltype)) +
  geom_point(aes(color = datatype))





