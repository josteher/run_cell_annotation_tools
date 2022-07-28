library(scID)
# Runs too slow also on smal datasets

## To do train labels as list


# Run scID 
# Input Data is raw counts
run_scID <- function(train, 
                     test, 
                     train_labels, 
                     test_labels){
  
  # Input formatting for method
  new_train_labels <- train_labels %>% set_names(colnames(train)) %>% as.factor()
  
  
  # Predict cell labels and measure time
  start_time <- Sys.time()
  
  scID_output <- scID::scid_multiclass(target_gem          = test, 
                                       reference_gem       = train, 
                                       reference_clusters  = new_train_labels,
                                       normalize_reference = TRUE)
  
  
  # Predictions of scID
  pred_labels <- scID_output$labels
  
  
  end_time <- Sys.time()
  runtime <- end_time - start_time
  
  
  return(list("output"     = scID_output,
              "runtime"    = runtime,
              "prediciton" = list("test_labels"  = test_labels,
                                  "pred_labels"  = pred_labels))
  )
}


# scID_test <- run_scID(train = TM_norm[,1:2000],
#                       test  = TM_norm[,2005:3000],
#                       train_labels = TM_labels[1:2000],
#                       test_labels  = TM_labels[2005:3000])

# scID_test <- run_scID(train = TM_norm[ , inds$train], 
#                       test  = TM_norm[ , inds$test], 
#                       train_labels = TM_labels[inds$train], 
#                       test_labels  = TM_labels[inds$test])



# # Example run
# scID_results <- run_scID(train = train, 
#                          test  = test, 
#                          train_labels = train_labels, 
#                          test_labels  = test_labels)
# 
# scID_results <- run_scID(train = train %>% as.matrix(), 
#                       test  = test %>% as.matrix(), 
#                       train_labels = train_meta$cell_ontology_class, 
#                       test_labels  = test_meta$cell_ontology_class)
# 



