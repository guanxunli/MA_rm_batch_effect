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

## simulation 2
#B cell (1199,1172)， CD4 T cell (2267,2183), CD8 T cell(2076, 1066), Monocyte_CD14 (1914, 2176)
cell_names = c("B cell", "CD4 T cell", "CD8 T cell", "Monocyte_CD14")
n1_sub = c(1000,1000,300,300)
n2_sub = c(300,300,1000,1000)

index1 = NA
index2 = NA
group_label1  = NA
group_label2  = NA
n1 = dim(b1_celltype)[1]
n2 = dim(b2_celltype)[1]
set.seed(123)
for(i in 1:4){
  index1 = c(index1, sample(c(1:n1)[c(b1_celltype$cell_type %in% cell_names[i])],n1_sub[i]))
  index2 = c(index2, sample(c(1:n2)[c(b2_celltype$cell_type %in% cell_names[i])],n2_sub[i]))
  group_label1 = c(group_label1, rep(i, n1_sub[i]))
  group_label2 = c(group_label2, rep(i, n2_sub[i]))
}

index1 = index1[-1]
index2 = index2[-1]
group_label1 = group_label1[-1]
group_label2 = group_label2[-1]

b1_exprs_sim2 = b1_exprs[, index1]
b2_exprs_sim2 = b2_exprs[, index2]
expr_mat_sim2 = cbind(b1_exprs_sim2,b2_exprs_sim2)
metadata_sim2 = rbind(b1_metadata[index1,], b2_metadata[index2,])
metadata_sim2$cell_type_index = c(group_label1, group_label2)



##########################################################
# run pipeline

source("R_Code/harmony/harmony_utility.R")
expr_mat = expr_mat_sim1
metadata = metadata_sim1
# data preprocess
npcs = 100
batch_label = "batchlb"
celltype_label = "cell_type"

b_seurat = data_preprocess(obj = expr_mat, min.cells = 10, min.features = 300, percent.mt = 10, oversd = NULL, 
                           normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst", 
                           nfeatures = 2000, npcs = npcs, metadata = metadata)

# plot before removing batch effect
b_seurat <- RunTSNE(b_seurat, reduction = "pca", seed.use = 10, dim.embed = 2, dims = 1:100)
b_seurat <- RunUMAP(b_seurat, reduction = "pca", n.components = 2, seed.use = 10 , dims = 1:100)
DimPlot(object = b_seurat, dims = c(1,2), reduction = "tsne", pt.size = 0.5, group.by = batch_label)
DimPlot(object = b_seurat, dims = c(1,2), reduction = "umap", pt.size = 0.5, group.by = batch_label)


##########################################################
# harmony in Seurat                           
b_seurat = RunHarmony(b_seurat, group.by.vars = batch_label,
                      theta_harmony = 2, numclust = 50, max_iter_cluster = 100)

b_seurat <- RunTSNE(b_seurat, reduction = "harmony", seed.use = 10, dim.embed = 2, dims = 1:100)
b_seurat <- RunUMAP(b_seurat, reduction = "harmony", n.components = 2, seed.use = 10 , dims = 1:100)
p11 <- DimPlot(object = b_seurat, dims = c(1,2), reduction = "tsne", pt.size = 0.5, group.by = batch_label)
p12 <- DimPlot(object = b_seurat, dims = c(1,2), reduction = "tsne", pt.size = 0.5, group.by = celltype_label)
print(p11 + p12)
p21 <- DimPlot(object = b_seurat, dims = c(1,2), reduction = "umap", pt.size = 0.5, group.by = batch_label)
p22 <- DimPlot(object = b_seurat, dims = c(1,2), reduction = "umap", pt.size = 0.5, group.by = celltype_label)
print(p21 + p22)

# save results
png("R_Code/results/harmony_results/pbmc/pbmc_tsne_harmony_sim1.png",width = 2*1000, height = 800, res = 2*72)
print(p11 + p12)
dev.off()

png("R_Code/results/harmony_results/pbmc/pbmc_umap_harmony_sim1.png",width = 2*1000, height = 800, res = 2*72)
print(p21 + p22)
dev.off()

batch_harmony = b_seurat$batch
cell_type_harmony = b_seurat$cell_type
saveRDS(b_seurat@reductions$harmony@cell.embeddings, "R_Code/results/harmony_results/pbmc/pbmc_harmony_sim1.rds")
save(batch_harmony, cell_type_harmony, file = "R_Code/results/harmony_results/pbmc/harmony_metadata_sim1.Rdata")


