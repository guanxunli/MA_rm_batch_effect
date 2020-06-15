library(harmony)
library(cowplot)
library(Seurat)
library(magrittr)
library(ggplot2)

rm(list=ls())

##########################################################
# read data 

b1_exprs <- readRDS("data set/PBMC_dataset/data1.rds")
b2_exprs <- readRDS("data set/PBMC_dataset/data2.rds")
b1_celltype <- readRDS("data set/PBMC_dataset/data1_celltype.rds")
b2_celltype <- readRDS("data set/PBMC_dataset/data2_celltype.rds")
rownames(b1_celltype) <- gsub("-", ".", x = b1_celltype[, 1])
b1_celltype <- b1_celltype[, 2, drop = FALSE]
rownames(b2_celltype) <- gsub("-", ".", x = b2_celltype[, 1])
b2_celltype <- b2_celltype[, 2, drop = FALSE]

b1_celltype <- b1_celltype[colnames(b1_exprs), , drop = FALSE]
b2_celltype <- b2_celltype[colnames(b2_exprs), , drop = FALSE]
b1_metadata <- as.data.frame(b1_celltype)
b2_metadata <- as.data.frame(b2_celltype)
b1_metadata$cell <- rownames(b1_celltype)
b2_metadata$cell <- rownames(b2_celltype)
b1_metadata$batch <- 1
b2_metadata$batch <- 2
b1_metadata$batchlb <- 'Batch_1'
b2_metadata$batchlb <- 'Batch_2'

expr_mat = cbind(b1_exprs,b2_exprs)
metadata = rbind(b1_metadata, b2_metadata)

expr_mat <- expr_mat[, rownames(metadata)]


#####################################################
## simulation 1 remove "NK cell" and "Monocyte_FCGR3A" in two dataset respectively
b1_exprs_sim1 = b1_exprs[, b1_celltype != "NK cell"]
b2_exprs_sim1 = b2_exprs[, b2_celltype != "Monocyte_FCGR3A"]
expr_mat_sim1 = cbind(b1_exprs_sim1,b2_exprs_sim1)
metadata_sim1 = rbind(b1_metadata[b1_celltype != "NK cell",], b2_metadata[b2_celltype != "Monocyte_FCGR3A",])

## run seurat
##########################################################
# run pipeline

source("R_Code/Seurat/Seurat_utility.R")

# data preprocess
npcs = 100
batch_label = "batchlb"
celltype_label = "cell_type"

batch_list <- data_preprocess(expr_mat_sim1, min.cells = 10, min.features = 300, percent.mt = 10, oversd = NULL, 
                              normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst", 
                              nfeatures = 2000, npcs = 50, metadata = metadata_sim1, batch_label = batch_label)

##########################################################
# Integration in Seurat  

batches <- seurat_integrate(batch_list = batch_list, npcs = npcs)

## plot results
dta_use <- as.matrix(batches@assays$integrated@data)
dta_use <- t(dta_use)
library(irlba)
dta_use_svd <- irlba(dta_use, nv = 100)
dta_use <- dta_use %*% dta_use_svd$v

tsne_embeddings <- Rtsne::Rtsne(dta_use, is_distance=FALSE, perplexity=30, num_threads=1, verbose=FALSE)$Y
umap_embeddings <- umap::umap(dta_use)$layout

p11 <- ggplot(data = NULL, aes(tsne_embeddings[, 1], tsne_embeddings[, 2], color = batches@meta.data$batchlb)) +
  geom_point(size = 0.5) +
  labs(x = "T-SNE1", y = "T-SNE2", color = "data type") 
p12 <- ggplot(data = NULL, aes(tsne_embeddings[, 1], tsne_embeddings[, 2], color = batches@meta.data$cell_type)) +
  geom_point(size = 0.5) +
  labs(x = "T-SNE1", y = "T-SNE2", color = "data type") 
print(plot_grid(p11, p12))
p21 <- ggplot(data = NULL, aes(umap_embeddings[, 1], umap_embeddings[, 2], color = batches@meta.data$batchlb)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP1", y = "UMAP2", color = "data type")
p22 <- ggplot(data = NULL, aes(umap_embeddings[, 1], umap_embeddings[, 2], color = batches@meta.data$cell_type)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP1", y = "UMAP2", color = "data type")
print(plot_grid(p21 + p22))

# save results
png("R_Code/Seurat/seurat_results/Simulation1/tsne_seurat.png",width = 2*1000, height = 800, res = 2*72)
print(plot_grid(p11, p12))
dev.off()

png("R_Code/Seurat/seurat_results/Simulation1/umap_seurat.png",width = 2*1000, height = 800, res = 2*72)
print(plot_grid(p21, p22))
dev.off()
