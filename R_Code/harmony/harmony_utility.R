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

harmony_seurat <- function(obj, batch_label = "batchlb",
                           theta_harmony = 2, numclust = 50, max_iter_cluster = 100){
  
  t1 = Sys.time()
  
  obj <- RunHarmony(object = obj, group.by.vars = batch_label, theta = theta_harmony, 
                    nclust = numclust, max.iter.cluster = max_iter_cluster)
  t2 = Sys.time()
  print(t2-t1)
  
  return(obj)
}
  
harmony_dir <- function(obj, metadata = NULL, batch_label = "batchlb", theta = 2, numclust = 50, max_iter_cluster = 100){
  dta <- obj@reductions$pca@cell.embeddings
  meta <- obj@meta.data
  harmony_embeddings <- HarmonyMatrix(Seurat_pca, meta_data = meta, vars_use = batch_label, do_pca=FALSE, 
                                      theta = theta, nclust = numclust, max.iter.cluster = max_iter_cluster)
  return(harmony_embeddings)
}

plot_res <- function(obj, reduc_method = "pca", seed.use = 10, dims = 1:npcs,
                     dataset = NULL, batch_label = "batchlb", celltype_label = "cell_type"){
  #############################
  #preparing plots
  
  obj <- RunTSNE(obj, reduction = reduc_method, seed.use = seed.use, dim.embed = 2, dims = dims)
  obj <- RunUMAP(obj, reduction = reduc_method, n.components = 2, seed.use = seed.use , dims = dims)
  
  #############################
  #tSNE plot
  p11 <- DimPlot(object = obj, dims = c(1,2), reduction = "tsne", pt.size = 0.5, group.by = batch_label)
  if (is.null(celltype_label) == FALSE){
    p12 <- DimPlot(object = obj, dims = c(1,2), reduction = "tsne", pt.size = 0.5, group.by = celltype_label)
    png(paste0("R_Code/harmony/results/", dataset, "_tsne_",reduc_method,".png"),width = 2*1000, height = 800, res = 2*72)
    print(plot_grid(p11, p12))
    dev.off()
  } else{
    png(paste0("R_Code/harmony/results/", dataset, "_tsne_",reduc_method,".png"),width = 2*1000, height = 800, res = 2*72)
    print(p11)
    dev.off()
  }

  
  #############################
  #UMAP plot
  p11 <- DimPlot(object = obj, dims = c(1,2), reduction = "umap", pt.size = 0.5, group.by = batch_label)
  if (is.null(celltype_label) == FALSE){
    p12 <- DimPlot(object = obj, dims = c(1,2), reduction = "umap", pt.size = 0.5, group.by = celltype_label)
    png(paste0("R_Code/harmony/results/", dataset, "_umap_",reduc_method,".png"),width = 2*1000, height = 800, res = 2*72)
    print(plot_grid(p11, p12))
    dev.off()
  } else{
    png(paste0("R_Code/harmony/results/", dataset, "_umap_",reduc_method,".png"),width = 2*1000, height = 800, res = 2*72)
    print(p11)
    dev.off()
  }
}


