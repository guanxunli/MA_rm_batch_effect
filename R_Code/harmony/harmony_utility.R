data_preprocess <- function(obj, min.cells = 10, min.features = 300, percent.mt = 10, oversd = NULL, 
                            normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst", 
                            nfeatures = 2000, npcs = 50, metadata = NULL) {
  
  obj <- Seurat::CreateSeuratObject(counts = obj, min.cells = min.cells, min.features = min.features, meta.data = metadata)
  
  # Do quality control
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, cells = which(obj[["percent.mt"]] < percent.mt))
  
  if (is.null(oversd) == FALSE){
    library_size <- obj$nCount_RNA
    sd_library_size <- sd(library_size)
    mean_library_size <- mean(library_size)
    thres_library_size <- mean_library_size + oversd * sd_library_size
    obj <- subset(obj, cells = which(obj[["nCount_RNA"]] < thres_library_size))
  }
  
  # Normalize object
  obj <- Seurat::NormalizeData(obj, normalization.method = normalization.method, scale.factor = scale.factor)
  
  # Identification of highly variable features (feature selection)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = selection.method, nfeatures = nfeatures)
  
  # Scaling the data
  all.genes <- rownames(obj)
  obj <- Seurat::ScaleData(obj, features = all.genes)
  
  # Perform linear dimensional reduction
  obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(object = obj), verbose = FALSE, npcs = npcs)
  
  # return Seurat object
  return(obj)
}

# plot for atacrna
plot_onetoone <- function(dta_embed, batch){
  dta_use <- data.frame("x" = dta_embed[, 1], "y" = dta_embed[, 2], "batch" = batch)
  n <- nrow(dta_embed)/2
  dta_use$group <- c(1:n, 1:n)
  p <- ggplot(data = dta_use, aes(x, y, color = batch, group = group)) +
    geom_point(size = 0.5) +
    geom_line(color = "black", size = 0.25) +
    labs(x = "embed1", y = "embed2", color = "data type")
  return(p)
}


