collaboratorlibrary(tidyverse)
library(symphony)
library(anndata)
library(Matrix)


source("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_scmap.R")
source("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_SingleR.R")
source("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Symphony.R")
source("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_singleCellNet.R")
source("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_scID.R")
source("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Seurat_refmap.R")

# Load Data and transform for methods

# Read data
tab_muris <- anndata::read_h5ad("/data/gruen/herman/scDatasets/Tabula_muris/Tabula_muris.h5ad")
anno <- read.csv("/data/gruen/herman/scDatasets/Tabula_muris/annotations_droplet.csv")

# Add cell_id
anno$barcode  <- sub("\\w*_","", anno$cell)
anno$cell_id  <- paste0(anno$barcode, '-1', '_', anno$tissue, '-', anno$channel)

# collaborator's splits
train_split <- read.csv("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/collaborator_splits_TM/train.split.csv", header = F) %>% .[,"V1"] %>% anno[.,]
test_split <- read.csv("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/collaborator_splits_TM/test.split.csv", header = F) %>% .[,"V1"] %>% anno[.,]
valid_split <- read.csv("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/collaborator_splits_TM/val.split.csv", header = F) %>% .[,"V1"] %>% anno[.,]

# Specify hidden classes as Unlabelled
hidden_classes <- c( "B cell", "hepatocyte", "keratinocyte", "mesenchymal stem cell" )

# Use 50% as unlabelled, specify top 50% of cells by class
unlab_top_50percent <- train_split %>% 
                          group_by(cell_ontology_class) %>%
                          arrange(cell_ontology_class) %>%
                          slice_head(prop = 0.5) %>%
                          ungroup() %>% select(cell_id) %>%
                          purrr::flatten_chr()


train <- tab_muris$X[train_split$cell_id,] %>% Matrix::t() %>% as.matrix()
test  <- tab_muris$X[test_split$cell_id,] %>% Matrix::t() %>% as.matrix()
valid <- tab_muris$X[valid_split$cell_id,] %>% Matrix::t() %>% as.matrix()






cat("\n\n\n\n\n\nData loaded\n\n\n\n\n\n\n")
### 
### RUN METHODS
###
### Run cell mapping methods
data_sets <- c("TM_hidden_plus50p_unlabelled")

for (data_set in data_sets){
  
 
  train_labels <- train_split$cell_ontology_class %>% purrr::set_names(train_split$cell_id)
  train_labels[train_labels %in% c("", hidden_classes )] <- "Unlabelled"
  
  test_labels  <- test_split$cell_ontology_class %>% purrr::set_names(test_split$cell_id)
  test_labels[test_labels %in% c("", hidden_classes )] <- "Unlabelled"
  
  valid_labels <- valid_split$cell_ontology_class %>% purrr::set_names(valid_split$cell_id)
  valid_labels[valid_labels %in% c("", hidden_classes )] <- "Unlabelled"
  
  
  if ( data_set == "TM_hidden_plus50p_unlabelled"){
    
    train_labels <- train_split$cell_ontology_class %>% purrr::set_names(train_split$cell_id)
    train_labels[train_labels %in% c("", hidden_classes )] <- "Unlabelled"
    # as.data.frame(train_labels) %>% tibble::rownames_to_column("cell_id") %>% write.table(., file="Trainset_collaborator_splits_hidden_classes.csv", quote = F, sep = ",",row.names = F)
    train_labels[names(train_labels) %in% unlab_top_50percent] <- "Unlabelled"
    # as.data.frame(train_labels) %>% tibble::rownames_to_column("cell_id") %>% write.table(., file="Trainset_collaborator_splits_hidden_classes_and_50p_unlabelled.csv", quote = F, sep = ",", row.names = F)
    
    test_labels  <- test_split$cell_ontology_class %>% purrr::set_names(test_split$cell_id)
    test_labels[test_labels %in% c("", hidden_classes )] <- "Unlabelled"
    
    valid_labels <- valid_split$cell_ontology_class %>% purrr::set_names(valid_split$cell_id)
    valid_labels[valid_labels %in% c("", hidden_classes )] <- "Unlabelled"
    
    
  }

  ##scmapcell
  scmapcell_results <- run_scmap(train, 
                                 test, 
                                 train_labels, 
                                 test_labels)
  saveRDS(scmapcell_results, file = paste0("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/",data_set,"_scmapcell_results.rds"))
  
  
  
  ## Symphony 
  symphony_results <- run_Symphony(train, 
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
                                   vargenes_groups = "label",
                                   do_normalize = T,
                                   do_umap = T)
  saveRDS(symphony_results, file = paste0("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/",data_set,"_symphony_results.rds"))
  
  
  
  ## SingleR
  SingleR_results <- run_SingleR(train, 
                                 test, 
                                 train_labels, 
                                 test_labels, 
                                 genes="de")
  saveRDS(SingleR_results, file = paste0("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/",data_set,"_SingleR_results.rds"))
  
  
  
  ## singleCellNet
  singleCellNet_results <- run_singleCellNet(train,
                                             test, 
                                             train_labels, 
                                             test_labels)
  saveRDS(singleCellNet_results, file = paste0("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/",data_set,"_singleCellNet_results.rds"))
  
  ## Seurat reference mapping
  Seurat_results <- run_Seurat_refmap(train = train,
                                      test  = test,
                                      train_labels = train_labels,
                                      test_labels  = test_labels)
  saveRDS(Seurat_results, file = paste0("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/",data_set,"_Seurat_refmap_results.rds"))
  
  
  ##scID
  scID_results <- run_scID(train = train,
                           test  = test,
                           train_labels = train_labels,
                           test_labels  = test_labels)
  saveRDS(scID_results, file = paste0("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/",data_set,"_scID_results.rds"))

}
