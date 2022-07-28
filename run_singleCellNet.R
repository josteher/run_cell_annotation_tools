library(singleCellNet)


run_singleCellNet <- function(train, 
                              test, 
                              train_labels, 
                              test_labels){
  
  # Create metadata dataframe for method
  train_meta_df <- data.frame("cell"  = colnames(train),
                              "label" = train_labels,
                              row.names = colnames(train) ) 
  
  test_meta_df <- data.frame("cell"  = colnames(test),
                             "label" = test_labels,
                             row.names = colnames(test)
  )
  
  # Subset only common genes
  common_genes <- intersect(rownames(train), rownames(test))
  train <- train[common_genes,]
  test  <- test[common_genes,]
  
  # Predict cell labels and measure time
  start_time <- Sys.time()
  
  # Train classifier
  class_info <- scn_train(stTrain = train_meta_df, 
                          expTrain = train, 
                          nTopGenes = 10, 
                          nRand = 70, 
                          nTrees = 1000, 
                          nTopGenePairs = 25, 
                          dLevel = "label", 
                          colName_samp = "cell",
                          )
    
  # Predict query dataset
  class_score <- scn_predict(class_info[['cnProc']], 
                        test, 
                        nrand = 0)
  
  # classify query cells
  query_assign <- assign_cate(classRes = class_score, sampTab = test_meta_df, cThresh = 0.5) 
  
  
  end_time <- Sys.time()
  runtime <- end_time - start_time
  
  
  return(list("runtime"    = runtime,
              "class_info" = class_info,
              "query_class_score" = class_score,
              "query_assign_label" = query_assign,
              "prediciton" = list("test_labels" = test_labels,
                                  "pred_labels" = query_assign$category))
  )
  
}
























