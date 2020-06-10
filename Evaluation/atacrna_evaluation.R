library(Seurat)
## load data
dta_harmony <- readRDS("R_Code/results/harmony_results/atacrna/harmony_atacrna.rds")
dta_harmonyseurat <- readRDS("R_Code/results/harmony_results/atacrna/harmonyseurat_atacrna.rds")
dta_seurat <- readRDS("R_Code/results/seurat_results/atacrna/seurat_atacrna.rds")

batch <- c(rep("atac", 1047), rep("rna", 1047))
library(irlba)
library(Matrix)
dta_seurat_svd <- irlba(dta_seurat, nv = 100)
dta_seurat <- dta_seurat %*% dta_seurat_svd$v

## evaluation
source("Evaluation/utility.R")
# kBET
dta_harmony_kBET <- kBET_fun(data = dta_harmony, batch = batch)
dta_harmonyseurat_kBET <- kBET_fun(data = dta_harmonyseurat, batch = batch)
dta_seurat_kBET <- kBET_fun(data = dta_seurat, batch = batch)
dta <- data.frame("reject rate" = c(dta_harmony_kBET, dta_harmonyseurat_kBET, dta_seurat_kBET), 
                     "percentage of sample size" = rep(seq(0.05, 0.25, by = 0.05), 3),
                     "data type" = c(rep("harmony-seperate", 5), rep("harmony-together", 5), rep("seurat", 5)))
ggplot(data = dta, aes(x = percentage.of.sample.size, y = reject.rate)) +
  geom_line(aes(color = data.type)) +
  geom_point(aes(color = data.type))

