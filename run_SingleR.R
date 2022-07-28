library(SingleR)

# Run SingleR with or withou variabele genes
# Input Data should be log transformed (e.g. using Seurat)
run_SingleR <- function(train, 
                        test, 
                        train_labels, 
                        test_labels, 
                        genes="de"){

# Predict cell labels and measure time
start_time <- Sys.time()
pred_labels <-  SingleR(test = test, 
                       ref = train, 
                       labels = train_labels,
                       genes = genes)

end_time <- Sys.time()
runtime <- end_time - start_time


return(list("runtime"    = runtime,
            "prediciton" = list("test_labels" = test_labels,
                                "pred_labels" = pred_labels))
       )

}


# 
# # Run Crossvalidtion for SingleR
# # Using fold created with splitTools
# cv_SingleR <- function(data, labels, folds, k_folds = 5){
#   
#   cv <- list()
#   # Run Crossvalidation
#   for (fold in folds){
#     cv[[paste0("Fold_", fold)]] <- run_SingleR(train = data[ , fold], 
#                                                test = data[ , -fold], 
#                                                train_labels = labels[fold], 
#                                                test_labels  = labels[-fold], 
#                                                genes="de")
#     return(cv)
#   }
# }

# 
# 
# ### Example run
# SingleR_results <- run_SingleR(train,
#                                test,
#                                train_labels,
#                                test_labels,
#                                genes="de")
# 
# SingleR_results <- run_SingleR(train = train %>% as.matrix(), 
#                                test  = test %>% as.matrix(), 
#                                train_labels = train_meta$cell_ontology_class, 
#                                test_labels  = test_meta$cell_ontology_class, 
#                                genes="de")











