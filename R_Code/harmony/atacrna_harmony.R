library(harmony)
library(cowplot)
library(Seurat)
library(magrittr)
library(ggplot2)

rm(list=ls())

##########################################################
# read data 

dta_atac <- read.csv("data set/dna_atac/dta_atac.csv")
atac_name <- dta_atac$X
dta_atac <- dta_atac[, -1]
rownames(dta_atac) <- atac_name

dta_rna <- read.csv("data set/dna_atac/dta_rna.csv")
rna_name <- dta_rna$X
dta_rna <- dta_rna[, -1]
rownames(dta_rna) <- rna_name

atac_metadata <- data.frame("batchlb" = rep("atac", ncol(dta_atac)))
rna_metadata <- data.frame("batchlb" = rep("rna", ncol(dta_rna)))

expr_mat = cbind(dta_atac, dta_rna)
metadata = rbind(atac_metadata, rna_metadata)
rownames(metadata) <- colnames(expr_mat)

##########################################################
# run pipeline

source("R_Code/harmony/harmony_utility.R")

# data preprocess
npcs = 100
batch_label = "batchlb"

b_seurat = data_preprocess(obj = expr_mat, min.cells = 1, min.features = 10, percent.mt = 10, oversd = NULL, 
                           normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst", 
                           nfeatures = 2000, npcs = npcs, metadata = metadata)

# plot before removing batch effect
b_seurat <- RunTSNE(b_seurat, reduction = "pca", seed.use = 10, dim.embed = 2, dims = 1:100)
b_seurat <- RunUMAP(b_seurat, reduction = "pca", n.components = 2, seed.use = 10 , dims = 1:100)
DimPlot(object = b_seurat, dims = c(1,2), reduction = "tsne", pt.size = 0.5, group.by = batch_label)
DimPlot(object = b_seurat, dims = c(1,2), reduction = "umap", pt.size = 0.5, group.by = batch_label)

# # save results
# plot_res(obj = b_seurat, reduc_method = "pca", dataset = "atac_rna", celltype_label = NULL)

##########################################################
# harmony in Seurat                           
b_seurat = harmony_seurat(b_seurat, batch_label = batch_label,
                          theta_harmony = 2, numclust = 50, max_iter_cluster = 100)
# plot before removing batch effect
b_seurat <- RunTSNE(b_seurat, reduction = "harmony", seed.use = 10, dim.embed = 2, dims = 1:100)
b_seurat <- RunUMAP(b_seurat, reduction = "harmony", n.components = 2, seed.use = 10 , dims = 1:100)
DimPlot(object = b_seurat, dims = c(1,2), reduction = "tsne", pt.size = 0.5, group.by = batch_label)
DimPlot(object = b_seurat, dims = c(1,2), reduction = "umap", pt.size = 0.5, group.by = batch_label)

# save results
plot_res(obj = b_seurat, reduc_method = "harmony", dataset = "atac_rna", celltype_label = NULL)
saveRDS(b_seurat@reductions$harmony@cell.embeddings, "R_Code/results/harmony_results/atacrna/harmonyseurat_atacrna.rds")
##################### seperate harmony ###################
rownames(rna_metadata) <- colnames(dta_rna)
b_seurat_rna <- data_preprocess(obj = dta_rna, min.cells = 1, min.features = 10, percent.mt = 10, oversd = NULL, 
                                normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst", 
                                nfeatures = 2000, npcs = npcs, metadata = rna_metadata)
rownames(atac_metadata) <- colnames(dta_atac)
b_seurat_atac <- data_preprocess(obj = dta_atac, min.cells = 1, min.features = 10, percent.mt = 10, oversd = NULL, 
                                normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst", 
                                nfeatures = 2000, npcs = npcs, metadata = atac_metadata)
dta_use <- rbind(b_seurat_atac@reductions$pca@cell.embeddings, b_seurat_rna@reductions$pca@cell.embeddings)
tsne_embeddings <- Rtsne::Rtsne(dta_use, is_distance=FALSE, perplexity=30, num_threads=1, verbose=FALSE)$Y
umap_embeddings <- umap::umap(dta_use)$layout
ggplot(data = NULL, aes(tsne_embeddings[, 1], tsne_embeddings[, 2], color = b_seurat@meta.data$batchlb)) +
  geom_point(size = 0.5) +
  labs(x = "T-SNE1", y = "T-SNE2", color = "data type")
ggplot(data = NULL, aes(umap_embeddings[, 1], umap_embeddings[, 2], color = b_seurat@meta.data$batchlb)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP1", y = "UMAP2", color = "data type")
## plot results
harmony_embeddings <- HarmonyMatrix(dta_use, meta_data = b_seurat@meta.data$batchlb, do_pca=FALSE)
tsne_embeddings <- Rtsne::Rtsne(harmony_embeddings, is_distance=FALSE, perplexity=30, num_threads=1, verbose=FALSE)$Y
umap_embeddings <- umap::umap(harmony_embeddings)$layout
ggplot(data = NULL, aes(tsne_embeddings[, 1], tsne_embeddings[, 2], color = b_seurat@meta.data$batchlb)) +
  geom_point(size = 0.5) +
  labs(x = "T-SNE1", y = "T-SNE2", color = "data type") 
ggplot(data = NULL, aes(umap_embeddings[, 1], umap_embeddings[, 2], color = b_seurat@meta.data$batchlb)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP1", y = "UMAP2", color = "data type")
## save results
saveRDS(harmony_embeddings, "R_Code/results/harmony_results/atacrna/harmony_atacrna.rds")
