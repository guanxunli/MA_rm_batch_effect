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

## run harmony
##########################################################
# run pipeline

source("R_Code/harmony/harmony_utility.R")

# data preprocess
npcs = 100
batch_label = "batchlb"
celltype_label = "cell_type"

b_seurat = data_preprocess(obj = expr_mat_sim1, min.cells = 10, min.features = 300, percent.mt = 10, oversd = NULL, 
                           normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst", 
                           nfeatures = 2000, npcs = npcs, metadata = metadata_sim1)

# plot before removing batch effect
b_seurat <- RunTSNE(b_seurat, reduction = "pca", seed.use = 10, dim.embed = 2, dims = 1:100)
b_seurat <- RunUMAP(b_seurat, reduction = "pca", n.components = 2, seed.use = 10 , dims = 1:100)
p1 <- DimPlot(object = b_seurat, dims = c(1,2), reduction = "tsne", pt.size = 0.5, group.by = batch_label)
p2 <- DimPlot(object = b_seurat, dims = c(1,2), reduction = "tsne", pt.size = 0.5, group.by = celltype_label)
plot(p1 + p2)
png("R_Code/harmony/harmony_results/Simulation1/sim1_tsne_pca.png",width = 2*1000, height = 800, res = 2*72)
plot(p1 + p2)
dev.off()

p1 <- DimPlot(object = b_seurat, dims = c(1,2), reduction = "umap", pt.size = 0.5, group.by = batch_label)
p2 <- DimPlot(object = b_seurat, dims = c(1,2), reduction = "umap", pt.size = 0.5, group.by = celltype_label)
plot(p1 + p2)
png("R_Code/harmony/harmony_results/Simulation1/sim1_umap_pca.png",width = 2*1000, height = 800, res = 2*72)
plot(p1 + p2)
dev.off()

##########################################################
# harmony in Seurat                           
b_seurat = RunHarmony(b_seurat, group.by.vars = batch_label,
                      theta_harmony = 2, numclust = 50, max_iter_cluster = 100)

b_seurat <- RunTSNE(b_seurat, reduction = "harmony", seed.use = 10, dim.embed = 2, dims = 1:100)
b_seurat <- RunUMAP(b_seurat, reduction = "harmony", n.components = 2, seed.use = 10 , dims = 1:100)
p11 <- DimPlot(object = b_seurat, dims = c(1,2), reduction = "tsne", pt.size = 0.5, group.by = batch_label)
p12 <- DimPlot(object = b_seurat, dims = c(1,2), reduction = "tsne", pt.size = 0.5, group.by = celltype_label)
plot(p11 + p12)
p21 <- DimPlot(object = b_seurat, dims = c(1,2), reduction = "umap", pt.size = 0.5, group.by = batch_label)
p22 <- DimPlot(object = b_seurat, dims = c(1,2), reduction = "umap", pt.size = 0.5, group.by = celltype_label)
plot(p21 + p22)

# save results
png("R_Code/harmony/harmony_results/Simulation1/sim1_tsne_harmony.png",width = 2*1000, height = 800, res = 2*72)
plot(p11 + p12)
dev.off()

png("R_Code/harmony/harmony_results/Simulation1/sim1_umap_harmony.png",width = 2*1000, height = 800, res = 2*72)
plot(p21 + p22)
dev.off()
