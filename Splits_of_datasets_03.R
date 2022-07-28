library(splitTools)
library(tidyverse)
# # 
# In order to preserve the saved csv files
# ONLY rerun if really needed

setwd("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/data_splits_12-07-22")

# Function for splitting the data
split_data <- function(cell_id,
                       labels,
                       train = 0.64, test = 0.2, valid = 0.16,
                       seed = 42,
                       # fold
                       k = 5,
                       bool_as_int = F,
                       hidden_classes = NULL
){

  split_stratified <- splitTools::partition(labels,
                                            p = c(train = train, test = test, valid = valid),
                                            seed = seed,
                                            type = "stratified")


  folds_stratified <- splitTools::create_folds(labels,
                                               k = k,
                                               seed = seed,
                                               type = "stratified")
  # split using different train & test sizes
  diff_size_splits <- list()

  # Split into train and test before splitting again to have one test set that is equal for all parameters
  pre_split <- splitTools::partition(labels,
                                     p = c(train = 0.8, test = 0.2),
                                     seed = seed,
                                     type = "stratified")

  for (i in seq(0.1, 0.9, 0.1) ){



    split <- splitTools::partition(pre_split$train,
                                   p = c(train = i, test = 1 - i),
                                   seed = seed,
                                   type = "stratified")
    split$train <- pre_split$train[split$train]
    split$test  <- pre_split$train[split$test]

    diff_size_splits[[paste0( "train_size_", as.character(i))]] <- split
  }


  split_df <- data.frame("cell_id" = cell_id,
                         "label" = labels,

                         "train_stratified" = seq_along(labels) %in% split_stratified$train,
                         "test_stratified"  = seq_along(labels) %in% split_stratified$test,
                         "valid_stratified" = seq_along(labels) %in% split_stratified$valid,

                         "Fold1" = seq_along(labels) %in% folds_stratified$Fold1,
                         "Fold2" = seq_along(labels) %in% folds_stratified$Fold2,
                         "Fold3" = seq_along(labels) %in% folds_stratified$Fold3,
                         "Fold4" = seq_along(labels) %in% folds_stratified$Fold4,
                         "Fold5" = seq_along(labels) %in% folds_stratified$Fold5,

                         "train_0.1" = seq_along(labels) %in% diff_size_splits$train_size_0.1$train,
                         "train_0.2" = seq_along(labels) %in% diff_size_splits$train_size_0.2$train,
                         "train_0.3" = seq_along(labels) %in% diff_size_splits$train_size_0.3$train,
                         "train_0.4" = seq_along(labels) %in% diff_size_splits$train_size_0.4$train,
                         "train_0.5" = seq_along(labels) %in% diff_size_splits$train_size_0.5$train,
                         "train_0.6" = seq_along(labels) %in% diff_size_splits$train_size_0.6$train,
                         "train_0.7" = seq_along(labels) %in% diff_size_splits$train_size_0.7$train,
                         "train_0.8" = seq_along(labels) %in% diff_size_splits$train_size_0.8$train,
                         "train_0.9" = seq_along(labels) %in% diff_size_splits$train_size_0.9$train,
                         "fixed_test_0.2" = seq_along(labels) %in% pre_split$test
  )
  if ( bool_as_int ){
    split_df <- split_df %>% mutate(across("train_stratified" : "fixed_test_0.2", as.integer))
  }

  if ( !is.null(hidden_classes) ){
    split_df <- split_df %>%
      mutate(new_label = ifelse(label %in% hidden_classes, "Unlabelled", label) )
  }

  return(split_df)
}





###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
### Tabula muris
###
TM_anno <- read.csv("/data/gruen/herman/scDatasets/Tabula_muris/annotations_droplet.csv")

# Add cell_id
TM_anno$barcode  <- sub("\\w*_","", TM_anno$cell)
TM_anno$cell_id  <- paste0(TM_anno$barcode, '-1', '_', TM_anno$tissue, '-', TM_anno$channel)

my_hidden_classes_TM <- c( "B cell", "hepatocyte", "keratinocyte", "mesenchymal stem cell", "" )


TM_split <- split_data(TM_anno$cell_id,
                       TM_anno$cell_ontology_class,
                       hidden_classes = my_hidden_classes_TM,
                       bool_as_int = T)


unlab_top_50percent <- TM_split %>%
                          group_by(label) %>%
                          arrange(label) %>%
                          slice_head(prop = 0.5) %>%
                          ungroup() %>% select(cell_id) %>%
                          purrr::flatten_chr()

TM_split <- TM_split %>%
              mutate(new_label_50p_unlab = ifelse(cell_id %in% unlab_top_50percent, "Unlabelled", new_label) )


write.table(TM_split, "Tabula_muris_split.csv", quote = F, sep = ",", row.names = T, col.names = T)

###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
### Tabula muris senis
###
TMS_anno <- read.csv("/data/gruen/herman/scDatasets/Tabula_muris_senis/tabula_muris_senis.obs.csv")
TMS_split <- split_data(TMS_anno$cell,
                        TMS_anno$cell_ontology_class,
                        bool_as_int = T)

write.table(TMS_split, "Tabula_muris_senis_split.csv", quote = F, sep = ",", row.names = T, col.names = T)

###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
### Tabula sapiens
###
TSap_anno <- read.csv("/data/gruen/herman/scDatasets/Tabula_sapiens/Tabula_sapiens_obs-annotation.csv")

TSap_split <- split_data(TSap_anno$cell_id,
                         TSap_anno$cell_ontology_class,
                         bool_as_int = T)

write.table(TSap_split, "Tabula_sapiens_split.csv", quote = F, sep = ",", row.names = T, col.names = T)


###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
### PBMC
###

pbmc_labels_10Xv2 <- readr::read_csv("/data/gruen/herman/scDatasets/Abdelaal_et_al_2019_Benchmarking/Inter-dataset/PbmcBench/10Xv2/10Xv2_pbmc1Labels.csv")$x
pbmc_labels_10Xv3 <- readr::read_csv("/data/gruen/herman/scDatasets/Abdelaal_et_al_2019_Benchmarking/Inter-dataset/PbmcBench/10Xv3/10Xv3_pbmc1Labels.csv")$x
pbmc_labels_celsq <- readr::read_csv("/data/gruen/herman/scDatasets/Abdelaal_et_al_2019_Benchmarking/Inter-dataset/PbmcBench/CEL-Seq/CL_pbmc1Labels.csv")$x
pbmc_labels_drops <- readr::read_csv("/data/gruen/herman/scDatasets/Abdelaal_et_al_2019_Benchmarking/Inter-dataset/PbmcBench/Drop-Seq/DR_pbmc1Labels.csv")$x
pbmc_labels_indrp <- readr::read_csv("/data/gruen/herman/scDatasets/Abdelaal_et_al_2019_Benchmarking/Inter-dataset/PbmcBench/inDrop/iD_pbmc1Labels.csv")$x
pbmc_labels_sqwel <- readr::read_csv("/data/gruen/herman/scDatasets/Abdelaal_et_al_2019_Benchmarking/Inter-dataset/PbmcBench/Seq-Well/SW_pbmc1Labels.csv")$x
pbmc_labels_smsq2 <- readr::read_csv("/data/gruen/herman/scDatasets/Abdelaal_et_al_2019_Benchmarking/Inter-dataset/PbmcBench/Smart-Seq2/SM2_pbmc1Labels.csv")$x


###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
### Pancreas
###



pc_baron_mouse_lab <- readr::read_csv("/data/gruen/herman/scDatasets/Abdelaal_et_al_2019_Benchmarking/Intra-dataset/Pancreatic_data/Baron Mouse/Labels.csv")$x
pc_baron_human_lab <- readr::read_csv("/data/gruen/herman/scDatasets/Abdelaal_et_al_2019_Benchmarking/Intra-dataset/Pancreatic_data/Baron Human/Labels.csv")$x
pc_muraro_lab <- readr::read_csv("/data/gruen/herman/scDatasets/Abdelaal_et_al_2019_Benchmarking/Intra-dataset/Pancreatic_data/Muraro/Labels.csv")$x
pc_segerstolp_lab <- readr::read_csv("/data/gruen/herman/scDatasets/Abdelaal_et_al_2019_Benchmarking/Intra-dataset/Pancreatic_data/Segerstolpe/Labels.csv")$x
pc_xin_lab <- data.table::fread("/data/gruen/herman/scDatasets/Abdelaal_et_al_2019_Benchmarking/Intra-dataset/Pancreatic_data/Xin/Labels.csv")$x



labels_list <- list("pbmc_labels_10Xv2" = pbmc_labels_10Xv2,
                    "pbmc_labels_10Xv3" = pbmc_labels_10Xv3,
                    "pbmc_labels_celseq" = pbmc_labels_celsq,
                    "pbmc_labels_dropseq" = pbmc_labels_drops,
                    "pbmc_labels_indrop" = pbmc_labels_indrp,
                    "pbmc_labels_seqwell" = pbmc_labels_sqwel,
                    "pbmc_labels_smartseq2" = pbmc_labels_smsq2,

                    "pc_baron_mouse_lab" = pc_baron_mouse_lab,
                    "pc_baron_human_lab" = pc_baron_human_lab,
                    "pc_muraro_lab" = pc_muraro_lab,
                    "pc_segerstolp_lab" = pc_segerstolp_lab,
                    "pc_xin_lab" = pc_xin_lab
                    )

labels_list_split <- lapply(labels_list, function(x){split_data(seq_along(x), x, bool_as_int = T)})


for ( i in names(labels_list_split)){
  write.table(labels_list_split[[i]], paste0(i, "_", ".csv"), quote = F, sep = ",", row.names = T, col.names = T)
}

