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

batch_list <- data_preprocess(expr_mat, min.cells = 10, min.features = 300, percent.mt = 10, oversd = NULL, 
                              normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst", 
                              nfeatures = 2000, npcs = 50, metadata = metadata, batch_label = batch_label)

##########################################################
# Integration in Seurat  

batches <- seurat_integrate(batch_list = batch_list, npcs = npcs)

plot_res(obj = batches, dataset = "atacrna", batch_label = batch_label, celltype_label = NULL)
