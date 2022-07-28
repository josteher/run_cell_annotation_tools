library(Seurat)

run_Seurat_refmap <- function(train, 
                              test, 
                              train_labels, 
                              test_labels){
  
  # Measure time
  start_time <- Sys.time()
  
  # Create reference Seurat object
  train.ref <- CreateSeuratObject(train)
  
  train.ref <- NormalizeData(train.ref)
  train.ref <- FindVariableFeatures(train.ref, 
                                   selection.method = "vst", 
                                   nfeatures = 2000,
                                   verbose = FALSE)
  
  train.ref <- ScaleData(train.ref, 
                        verbose = FALSE)
  train.ref <- RunPCA(train.ref, 
                     npcs = 30, 
                     verbose = FALSE)
  train.ref <- RunUMAP(train.ref, 
                      reduction = "pca", 
                      dims = 1:30, 
                      verbose = FALSE)
  
  # Create query Seurat object
  test.query <- CreateSeuratObject(test)
  test.query <- NormalizeData(test.query)
  test.query <- FindVariableFeatures(test.query, 
                                    selection.method = "vst", 
                                    nfeatures = 2000,
                                    verbose = FALSE)
  
  test.query <- ScaleData(test.query, 
                         verbose = FALSE)
  test.query <- RunPCA(test.query, 
                      npcs = 30, 
                      verbose = FALSE)
  test.query <- RunUMAP(test.query, 
                       reduction = "pca", 
                       dims = 1:30, 
                       verbose = FALSE)
  
  
  # Transfer data and predict cell types
  train.anchors <- FindTransferAnchors(reference = train.ref, 
                                       query = test.query,
                                       dims = 1:30, 
                                       reference.reduction = "pca")
  
  predictions <- TransferData(anchorset = train.anchors, 
                              refdata = train_labels,
                              dims = 1:30)
  test.query <- AddMetaData(test.query, 
                            metadata = predictions)
  # all(rownames(test.query@meta.data) == colnames(test)) # == TRUE
  
  Seurat_output <- list('ref_data'   = train.ref,
                       'query_data' = test.query)
  
  end_time <- Sys.time()
  runtime <- end_time - start_time
  
  return(list("output"     = Seurat_output,
              "runtime"    = runtime,
              "prediciton" = list("test_labels"  = test_labels,
                                  "pred_labels"  = test.query@meta.data$predicted.id))
  )
}

