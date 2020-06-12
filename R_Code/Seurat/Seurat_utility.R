data_preprocess <- function(obj, min.cells = 10, min.features = 300, percent.mt = 10, oversd = NULL, 
                            normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst", 
                            nfeatures = 2000, npcs = 50, metadata = NULL, batch_label = batch_label){
  
  ##########################################################
  # preprocessing
  batches <- CreateSeuratObject(obj, meta.data = metadata, 
                                min.cells = min.cells, min.features = min.features, project = "Seurat3_benchmark")
  
  # Do quality control
  batches[["percent.mt"]] <- Seurat::PercentageFeatureSet(batches, pattern = "^MT-")
  batches <- subset(batches, cells = which(batches[["percent.mt"]] < percent.mt))
  
  if (is.null(oversd) == FALSE){
    library_size <- batches$nCount_RNA
    sd_library_size <- sd(library_size)
    mean_library_size <- mean(library_size)
    thres_library_size <- mean_library_size + oversd * sd_library_size
    batches <- subset(batches, cells = which(batches[["nCount_RNA"]] < thres_library_size))
  }
  
  batch_list <- SplitObject(batches, split.by = batch_label)
  for(i in 1:length(batch_list)) {
    # Normalize object
    batch_list[[i]] <- Seurat::NormalizeData(batch_list[[i]], normalization.method = normalization.method, scale.factor = scale.factor)
    
    # Identification of highly variable features (feature selection)
    batch_list[[i]] <- Seurat::FindVariableFeatures(batch_list[[i]], selection.method = selection.method, nfeatures = nfeatures)
  }
  return(batch_list)
}

seurat_integrate <- function(batch_list, npcs = 50){
  t1 = Sys.time()
  
  cell_anchors <- FindIntegrationAnchors(object.list = batch_list, dims = 1:npcs)
  batches <- IntegrateData(anchorset = cell_anchors, dims = 1:npcs)
  
  t2 = Sys.time()
  print(t2-t1)
  
  DefaultAssay(batches) <- "integrated"
  batches <- ScaleData(batches, verbose = FALSE)
  batches <- RunPCA(batches, npcs = npcs, verbose = FALSE)
  
  return(batches)
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


