library(SCINA)




run_SCINA <- function(train, 
                     test, 
                     train_labels, 
                     test_labels){
  
  # Predict cell labels and measure time
  start_time <- Sys.time()
  pred_labels <-  
    
    
  end_time <- Sys.time()
  runtime <- end_time - start_time
  
  
  return(list("runtime"    = runtime,
              "prediciton" = list("test_labels" = test_labels,
                                  "pred_labels" = pred_labels))
  )
  
}















