library(scmap)
library(SingleCellExperiment)


run_scmap <- function(train, 
                     test, 
                     train_labels, 
                     test_labels){
  
  ### Input data formatting for scmap
  # Training dataset
  sce_train <- SingleCellExperiment(assays  = list(normcounts = as.matrix(train)), 
                                    colData = data.frame(row.names = colnames(train),
                                                         "cell_type1"  = train_labels ) )
  
  logcounts(sce_train) <- log2(normcounts(sce_train) + 1)
  # use gene names as feature symbols
  rowData(sce_train)$feature_symbol <- rownames(sce_train)

  ###
  # Test dataset
  sce_test <- SingleCellExperiment(assays  = list(normcounts = as.matrix(test)), 
                                   colData = data.frame(row.names = colnames(test),
                                                        "cell_type1"  = test_labels ) )
  logcounts(sce_test) <- log2(normcounts(sce_test) + 1)
  # use gene names as feature symbols
  rowData(sce_test)$feature_symbol <- rownames(sce_test)
  
  # Select Features
  sce_train <- selectFeatures(sce_train, suppress_plot = FALSE)
  
  ### Predict cell labels and measure time scmap cluster
  start_time <- Sys.time()
  # Index groups
  sce_train <- indexCluster(sce_train)
  
  # Project testset cells
  scmapCluster_results <- scmapCluster(
    projection = sce_test, 
    index_list = list(
      yan = metadata(sce_train)$scmap_cluster_index
    )
  )
  
  # Predicted Labels
  pred_labels_cluster <- scmapCluster_results$combined_labs
  end_time <- Sys.time()
  runtime <- end_time - start_time
  

  
  ##
  ### Predict cell labels and measure time scmap cell
  # i.e. projection on individual cells)
  start_time <- Sys.time()
  # Index groups
  sce_train <- indexCell(sce_train)
  
  # Projection
  scmapCell_results <- scmapCell(
    sce_test, 
    list(
      yan = metadata(sce_train)$scmap_cell_index
    )
  )
  
  # Cluster annotation
  scmapCell_clusters <- scmapCell2Cluster(
    scmapCell_results, 
    list(
      as.character(colData(sce_train)$cell_type1)
    )
  )
  
  # Predicted labels
  pred_labels_cell <- scmapCell_clusters$combined_labs
  
  end_time <- Sys.time()
  runtime2 <- end_time - start_time
  
  return(list("output" = 
                list("scmapCluster_results"     = scmapCluster_results,
                     "scmapCell_results"        = scmapCell_results,
                     "scmapCell_clusters"       = scmapCell_clusters),
              "runtime_scmap_cluster"    = runtime,
              "runtime_scmap_cell"       = runtime2,
              "test_labels" = test_labels,
              "prediciton" = 
                list("pred_labels_cluster" = pred_labels_cluster,
                     "pred_labels_cell"    = pred_labels_cell
                     )
              )
  )
  
}



# ## Example run
# scmap_test <- run_scmap(train = TM_filt_mat[ , inds$train], 
#                        test  = TM_filt_mat[ , inds$test], 
#                        train_labels = TM_labels[inds$train], 
#                        test_labels  = TM_labels[inds$test])
# 
# scmap_test <- run_scmap(train = train %>% as.matrix(), 
#                         test  = test %>% as.matrix(), 
#                         train_labels = train_meta$cell_ontology_class, 
#                         test_labels  = test_meta$cell_ontology_class)






















