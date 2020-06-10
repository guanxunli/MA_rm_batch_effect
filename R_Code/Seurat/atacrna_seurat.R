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

source("R_Code/Seurat/Seurat_utility.R")

# data preprocess
npcs = 50
batch_label = "batchlb"

batch_list <- data_preprocess(expr_mat, min.cells = 1, min.features = 10, percent.mt = 10, oversd = NULL, 
                              normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst", 
                              nfeatures = 2000, npcs = 50, metadata = metadata, batch_label = batch_label)

##########################################################
# Integration in Seurat  

batches <- seurat_integrate(batch_list = batch_list, npcs = npcs)

## plot results

## plot results
dta_use <- as.matrix(batches@assays$integrated@data)
dta_use <- t(dta_use)
tsne_embeddings <- Rtsne::Rtsne(dta_use, is_distance=FALSE, perplexity=30, num_threads=1, verbose=FALSE)$Y
umap_embeddings <- umap::umap(dta_use)$layout

p1 <- ggplot(data = NULL, aes(tsne_embeddings[, 1], tsne_embeddings[, 2], color = metadata$batchlb)) +
  geom_point(size = 0.5) +
  labs(x = "T-SNE1", y = "T-SNE2", color = "data type") 
print(p1)
p2 <- ggplot(data = NULL, aes(umap_embeddings[, 1], umap_embeddings[, 2], color = metadata$batchlb)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP1", y = "UMAP2", color = "data type")
print(p2)

# save results
png("R_Code/results/seurat_results/atacrna/tsne_seurat.png",width = 2*1000, height = 800, res = 2*72)
print(p1)
dev.off()

png("R_Code/results/seurat_results/atacrna/umap_seurat.png",width = 2*1000, height = 800, res = 2*72)
print(p2)
dev.off()

saveRDS(dta_use, "R_Code/results/seurat_results/atacrna/seurat_atacrna.rds")
