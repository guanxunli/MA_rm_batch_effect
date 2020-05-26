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

plot_res <- function(obj, seed.use = 10, dims = 1:npcs,
                     dataset = NULL, batch_label = batch_label, celltype_label = celltype_label){
  #############################
  #preparing plots
  
  obj <- RunTSNE(obj, reduction = "pca", seed.use = seed.use, dim.embed = 2, dims = dims)
  obj <- RunUMAP(obj, reduction = "pca", n.components = 2, seed.use = seed.use , dims = dims)
  
  #############################
  #tSNE plot
  p11 <- DimPlot(object = obj, dims = c(1,2), reduction = "tsne", pt.size = 0.5, group.by = batch_label)
  if (is.null(celltype_label) == FALSE){
    p12 <- DimPlot(object = obj, dims = c(1,2), reduction = "tsne", pt.size = 0.5, group.by = celltype_label)
    png(paste0("R_Code/results/seurat_results/", dataset, "/tsne_seurat.png"),width = 2*1000, height = 800, res = 2*72)
    print(plot_grid(p11, p12))
    dev.off()
  } else{
    png(paste0("R_Code/results/seurat_results/", dataset, "/tsne_seurat.png"),width = 2*1000, height = 800, res = 2*72)
    print(p11)
    dev.off()
  }
  
  #############################
  #UMAP plot
  p11 <- DimPlot(object = obj, dims = c(1,2), reduction = "umap", pt.size = 0.5, group.by = batch_label)
  if (is.null(celltype_label) == FALSE){
    p12 <- DimPlot(object = obj, dims = c(1,2), reduction = "umap", pt.size = 0.5, group.by = celltype_label)
    png(paste0("R_Code/results/seurat_results/", dataset, "/umap_seurat.png"),width = 2*1000, height = 800, res = 2*72)
    print(plot_grid(p11, p12))
    dev.off()
  } else{
    png(paste0("R_Code/results/seurat_results/", dataset, "/umap_seurat.png"),width = 2*1000, height = 800, res = 2*72)
    print(p11)
    dev.off()
  }
}


