library(symphony)
library(Matrix)


run_Symphony <- function(train, 
                         test, 
                         train_labels, 
                         test_labels,
                         ref_metadata = data.frame("cell"  = colnames(train),
                                                   "label" = train_labels,
                                                   "batch" = "batch1",
                                                   row.names = colnames(train) ),
                         query_metadata = data.frame("cell"  = colnames(test),
                                                     "label" = test_labels,
                                                     "batch" = "batch1",
                                                     row.names = colnames(test) ),
                         vargenes_groups = "batch",
                         topn = 2000,
                         do_normalize = T,
                         do_umap = T
){
  
  
  # train <- train %>% Matrix::t() %>% as(., "dgCMatrix")
  # test  <- test  %>% Matrix::t() %>% as(., "dgCMatrix")
  
  train <- train  %>% as(., "dgCMatrix")
  test  <- test   %>% as(., "dgCMatrix")
  
  
  start_time <- Sys.time()
  
  # Build reference
  set.seed(0)
  reference = symphony::buildReference(
    train,
    ref_metadata,
    vars = c('batch'),         # variables to integrate over
    K = 100,                   # number of Harmony clusters
    verbose = TRUE,            # verbose output
    do_umap = do_umap,            # can set to FALSE if want to run umap separately later
    do_normalize = do_normalize,      # set to TRUE if input counts are not normalized yet
    vargenes_method = 'vst',   # method for variable gene selection ('vst' or 'mvp')
    vargenes_groups = vargenes_groups, # metadata column specifying groups for variable gene selection 
    topn = topn,               # number of variable genes to choose per group
    d = 20#,                    # number of PCs
    # save_uwot_path = './testing_uwot_model_2'
  )
  
  
  # Map query
  query = mapQuery(test,             # query gene expression (genes x cells)
                   query_metadata,        # query metadata (cells x attributes)
                   reference,             # Symphony reference object
                   do_normalize = do_normalize,  # perform log(CP10k+1) normalization on query
                   do_umap = do_umap)        # project query cells into reference UMAP
  
  query = knnPredict(query, reference, reference$meta_data$label, k = 5)
  
  # Measure runtime
  end_time <- Sys.time()
  runtime <- end_time - start_time
  
  return(list("reference" = reference,
               "query"    = query,
               "runtime"  = runtime,
              "prediction" = list("test_labels" = test_labels,
                                  "pred_labels" = query$meta_data$cell_type_pred_knn) )
         )
  
  
  }