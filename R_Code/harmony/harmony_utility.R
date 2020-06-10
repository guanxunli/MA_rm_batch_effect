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



