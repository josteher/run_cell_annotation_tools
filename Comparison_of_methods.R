symphony_res <- readRDS("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/TM_hidden_symphony_results.rds")
scmapcell_res <- readRDS("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/TM_hidden_scmapcell_results.rds")
singler_res <- readRDS("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/TM_hidden_SingleR_results.rds")
singleCellNet_res <- readRDS("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/TM_hidden_singleCellNet_results.rds")
Seurat_res <- readRDS("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/TM_hidden_Seurat_refmap_results.rds")
scID_res <- readRDS("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/TM_hidden_scID_results.rds")

symphony_res_50p <- readRDS("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/TM_hidden_plus50p_unlabelled_symphony_results.rds")
scmapcell_res_50p <- readRDS("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/TM_hidden_plus50p_unlabelled_scmapcell_results.rds")
singler_res_50p <- readRDS("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/TM_hidden_plus50p_unlabelled_SingleR_results.rds")
singleCellNet_res_50p <- readRDS("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/TM_hidden_plus50p_unlabelled_singleCellNet_results.rds")
Seurat_res_50p <- readRDS("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/TM_hidden_plus50p_unlabelled_Seurat_refmap_results.rds")
scID_res_50p <- readRDS("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/TM_hidden_plus50p_unlabelled_scID_results.rds")


# Calculate class frequency
class_frequency <- read.csv("/data/gruen/herman/scDatasets/Tabula_muris/annotations_droplet.csv") %>% 
                dplyr::select(cell_ontology_class) %>% table() %>% sort(decreasing = T)
names(class_frequency)[names(class_frequency) == ""] <- "unassigned"
names(class_frequency) <- names(class_frequency) %>% gsub("\\s","_", .) %>% tolower()


# Summary df of predictions
res_df <- data.frame("true_test_label" = scmapcell_res$test_labels,
                     "symphony_pred"   = symphony_res$query$meta_data$cell_type_pred_knn,
                     "scmapcell_pred"  = scmapcell_res$prediciton$pred_labels_cell,
                     "singler_pred"    = singler_res$prediciton$pred_labels$labels,
                     "singcenet_pred"  = singleCellNet_res$prediciton$pred_labels,
                     "seurat_pred"     = Seurat_res$prediciton$pred_labels,
                     "scid_pred"       = scID_res$prediciton$pred_labels,
                     
                     "symphony_50p_pred"   = symphony_res_50p$query$meta_data$cell_type_pred_knn,
                     "scmapcell_50p_pred"  = scmapcell_res_50p$prediciton$pred_labels_cell,
                     "singler_50p_pred"    = singler_res_50p$prediciton$pred_labels$labels,
                     "singcenet_50p_pred"  = singleCellNet_res_50p$prediciton$pred_labels,
                     "seurat_50p_pred"     = Seurat_res_50p$prediciton$pred_labels,
                     "scid_50p_pred"       = scID_res_50p$prediciton$pred_labels
                     )

write_csv(res_df, "/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/prediction_summary.csv")

# # Reformat classification report, which was saved as text

### NEEDS TO BE DONE once
# can be avoided:
# classification_reports should be exported as dict then converted to pandas data frame

# setwd("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/classification_reports_data")
# 
# writeLines(
#   c(("X,precision,recall,f1.score,support"),
#   read.csv("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/muris/fully_supervised/negative_binomial/entropy-b35164fb/classification_report.dat", 
#                                   sep = "\t") %>% .[-1,] %>% gsub("\\s", "_",.) %>% gsub("^_*", "",.) %>% gsub("[_]{2,20}", ",",.)),
#   "nb_fully_supervised_b35164fb.csv"
# )
# 
# 
# writeLines(
#   c(("X,precision,recall,f1.score,support"),
#     read.csv("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/muris/fully_supervised/negative_binomial/entropy-e2e24638/classification_report.dat", 
#              sep = "\t") %>% .[-1,] %>% gsub("\\s", "_",.) %>% gsub("^_*", "",.) %>% gsub("[_]{2,20}", ",",.)),
#   "nb_fully_supervised_e2e24638.csv"
# )
# 
# 
# writeLines(
#   c(("X,precision,recall,f1.score,support"),
#     read.csv("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/muris/fully_supervised/gauss/entropy-1b5d8d37/classification_report.dat", 
#              sep = "\t") %>% .[-1,] %>% gsub("\\s", "_",.) %>% gsub("^_*", "",.) %>% gsub("[_]{2,20}", ",",.)),
#   "gauss_fully_supervised_1b5d8d37.csv"
# )
# 
# 
# writeLines(
#   c(("X,precision,recall,f1.score,support"),
#     read.csv("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/muris/fully_supervised/gauss/entropy-f131e011/classification_report.dat", 
#              sep = "\t") %>% .[-1,] %>% gsub("\\s", "_",.) %>% gsub("^_*", "",.) %>% gsub("[_]{2,20}", ",",.)),
#   "gauss_fully_supervised_f131e011.csv"
# )
# 
# 
# 
# writeLines(
#   c(("X,precision,recall,f1.score,support"),
#     read.csv("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/muris/semi_supervised_50p/negative_binomial/entropy-11b0a4c5/classification_report.dat", 
#              sep = "\t") %>% .[-1,] %>% gsub("\\s", "_",.) %>% gsub("^_*", "",.) %>% gsub("[_]{2,20}", ",",.)),
#   "semi_supervised_50p_11b0a4c5.csv"
# )
# 
# 
# writeLines(
#   c(("X,precision,recall,f1.score,support"),
#     read.csv("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/muris/semi_unsupervised_50p/entropy-0ba64550/classification_report.dat", 
#              sep = "\t") %>% .[-1,] %>% gsub("\\s", "_",.) %>% gsub("^_*", "",.) %>% gsub("[_]{2,20}", ",",.)),
#   "semi_unsupervised_50p_0ba64550.csv"
# )
# 
# 
# writeLines(
#   c(("X,precision,recall,f1.score,support"),
#     read.csv("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/muris/semi_unsupervised_100p/entropy-12ffe350/classification_report.dat", 
#              sep = "\t") %>% .[-1,] %>% gsub("\\s", "_",.) %>% gsub("^_*", "",.) %>% gsub("[_]{2,20}", ",",.)),
#   "semi_unsupervised_100p_12ffe350.csv"
# )


##
### !!!! in the accuracy line two additional NA values need to be added manually
## classification_reports should be exported as dict then converted to pandas data frame





classification_report <- list("symphony" = read.csv(file = "/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/classification_reports/symphony_pred.csv",
                                                    header = T),
                              "scmapcell" = read.csv(file = "/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/classification_reports/scmapcell_pred.csv",
                                                    header = T),
                              "singler" = read.csv(file = "/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/classification_reports/singler_pred.csv",
                                                    header = T),
                              "singleCellNet" = read.csv(file = "/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/classification_reports/singcenet_pred.csv",
                                                                      header = T),
                              "Seurat"     = read.csv(file = "/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/classification_reports/seurat_pred.csv",
                                                      header = T),
                              "scID"       = read.csv(file = "/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/classification_reports/scid_pred.csv",
                                                      header = T),
                              
                              
                              "symphony_50p" = read.csv(file = "/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/classification_reports/symphony_50p_pred.csv",
                                                    header = T),
                              "scmapcell_50p" = read.csv(file = "/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/classification_reports/scmapcell_50p_pred.csv",
                                                     header = T),
                              "singler_50p" = read.csv(file = "/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/classification_reports/singler_50p_pred.csv",
                                                   header = T),
                              "singleCellNet_50p" = read.csv(file = "/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/classification_reports/singcenet_50p_pred.csv",
                                                         header = T),
                              "Seurat_50p"     = read.csv(file = "/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/classification_reports/seurat_50p_pred.csv",
                                                      header = T),
                              "scID_50p"       = read.csv(file = "/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/classification_reports/scid_50p_pred.csv",
                                                      header = T),
                              
                              
                              "X_nb_fully_supervised_1" = read.csv("nb_fully_supervised_b35164fb.csv", header = T),
                              "X_nb_fully_supervised_2" = read.csv("nb_fully_supervised_e2e24638.csv", header = T),
                              "X_gauss_fully_supervised_1" = read.csv("gauss_fully_supervised_1b5d8d37.csv", header = T),
                              "X_gauss_fully_supervised_2" = read.csv("gauss_fully_supervised_f131e011.csv", header = T),
                              "X_semi_supervised_50p" = read.csv("semi_supervised_50p_11b0a4c5.csv", header = T),
                              "X_semi_unsupervised_50p" = read.csv("semi_unsupervised_50p_0ba64550.csv", header = T),
                              "X_semi_unsupervised_100p" = read.csv("semi_unsupervised_100p_12ffe350.csv", header = T)
                              
                              ) %>% dplyr::bind_rows(.id = "method") 


classification_report$X <- gsub("_\\(hidden\\)", "", classification_report$X) %>% tolower() %>% gsub("\\s", "_",.)



class_report_long <- classification_report %>% 
                        pivot_longer(cols = precision:f1.score, names_to = "metric") %>% 
                        arrange(desc(value)) %>% 
                        filter(!X %in% c("accuracy", "macro_avg", "weighted_avg"))

class_report_sum <- classification_report %>% 
                        filter(X %in% c("accuracy", "macro_avg", "weighted_avg")) %>% 
                        pivot_longer(cols = precision:f1.score, names_to = "metric") 



# Precision
class_report_long %>% 
  filter(metric == "precision") %>% 
  filter(method %in% c("Seurat", "singleCellNet", "singler", "symphony","scmapcell","scID", "X_gauss_fully_supervised_1",
                       "X_nb_fully_supervised_1")) %>% 
  arrange(desc(value)) %>% 
  ggplot(aes( y = value, x = factor(X, levels=c(names(class_frequency),"unlabelled") ), color = method, label = round(value,3) ) ) +
  geom_point(alpha = 0.8, size = 2,) + 
  geom_text(size = 2,nudge_x = 0.5) +
  scale_colour_manual(values =  RColorBrewer::brewer.pal(8, "Paired")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 



class_report_long %>% 
  filter(metric == "precision") %>% 
  filter(method %in% c("Seurat_50p", "singleCellNet_50p", "singler_50p", "symphony_50p","scmapcell_50p","scID_50p", "X_semi_supervised_50p",
                       "X_semi_unsupervised_50p")) %>% 
  arrange(desc(value)) %>% 
  ggplot(aes( y = value, x = factor(X, levels=c(names(class_frequency),"unlabelled") ), color = method, label = round(value,3) ) ) +
  geom_point(alpha = 0.8, size = 2,) + 
  geom_text(size = 2, nudge_x = 0.5) +
  scale_colour_manual(values =  RColorBrewer::brewer.pal(8, "Paired")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 




# Recall
class_report_long %>% 
  filter(metric == "recall") %>% 
  filter(method %in% c("Seurat", "singleCellNet", "singler", "symphony","scmapcell","scID", "X_gauss_fully_supervised_1",
                       "X_nb_fully_supervised_1")) %>% 
  arrange(desc(value)) %>% 
  ggplot(aes( y = value, x = factor(X, levels=c(names(class_frequency),"unlabelled") ), color = method, label = round(value,3) ) ) +
  geom_point(alpha = 0.8, size = 2,) + 
  geom_text(size = 2,nudge_x = 0.5) +
  scale_colour_manual(values =  RColorBrewer::brewer.pal(8, "Paired")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 



class_report_long %>% 
  filter(metric == "recall") %>% 
  filter(method %in% c("Seurat_50p", "singleCellNet_50p", "singler_50p", "symphony_50p","scmapcell_50p","scID_50p", "X_semi_supervised_50p",
                       "X_semi_unsupervised_50p")) %>% 
  arrange(desc(value)) %>% 
  ggplot(aes( y = value, x = factor(X, levels=c(names(class_frequency),"unlabelled") ), color = method, label = round(value,3) ) ) +
  geom_point(alpha = 0.8, size = 2,) + 
  geom_text(size = 2,nudge_x = 0.5) +
  scale_colour_manual(values =  RColorBrewer::brewer.pal(8, "Paired")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 




# F1 score
class_report_long %>% 
  filter(metric == "f1.score") %>% 
  filter(method %in% c("Seurat", "singleCellNet", "singler", "symphony","scmapcell","scID", "X_gauss_fully_supervised_1",
                       "X_nb_fully_supervised_1")) %>% 
  arrange(desc(value)) %>% 
  ggplot(aes( y = value, x = factor(X, levels=c(names(class_frequency),"unlabelled") ), color = method, label = round(value,3) ) ) +
  geom_point(alpha = 0.8, size = 2,) + 
  geom_text(size = 2,nudge_x = 0.5) +
  scale_colour_manual(values =  RColorBrewer::brewer.pal(8, "Paired")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 



class_report_long %>% 
  filter(metric == "f1.score") %>% 
  filter(method %in% c("Seurat_50p", "singleCellNet_50p", "singler_50p", "symphony_50p","scmapcell_50p","scID_50p", "X_semi_supervised_50p",
                       "X_semi_unsupervised_50p")) %>% 
  arrange(desc(value)) %>% 
  ggplot(aes( y = value, x = factor(X, levels=c(names(class_frequency),"unlabelled") ), color = method, label = round(value,3) ) ) +
  geom_point(alpha = 0.8, size = 2,) + 
  geom_text(size = 2,nudge_x = 0.5) +
  scale_colour_manual(values =  RColorBrewer::brewer.pal(8, "Paired")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 







# Summary statistic

class_report_sum %>% 
  # filter(metric == "f1.score") %>% 
  filter(method %in% c("Seurat", "singleCellNet", "singler", "symphony","scmapcell","scID", "X_gauss_fully_supervised_1",
                       "X_nb_fully_supervised_1")) %>% 
  arrange(desc(value)) %>% 
  ggplot(aes( y = value, x = X, color = method, label = round(value,3) ) ) +
  geom_point(alpha = 0.8, size = 1.5,) + 
  ggrepel::geom_text_repel(size = 2.5) +
  scale_colour_manual(values =  RColorBrewer::brewer.pal(8, "Paired")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = -40, hjust = 0)) +
  facet_wrap(vars(metric))


class_report_sum %>% 
  # filter(metric == "f1.score") %>% 
  filter(method %in% c("Seurat_50p", "singleCellNet_50p", "singler_50p", "symphony_50p","scmapcell_50p","scID_50p", "X_semi_supervised_50p",
                       "X_semi_unsupervised_50p", "X_semi_unsupervised_100p")) %>% 
  arrange(desc(value)) %>% 
  ggplot(aes( y = value, x = X, color = method, label = round(value,3) ) ) +
  geom_point(alpha = 0.8, size = 1.5,) + 
  ggrepel::geom_text_repel(size = 2.5) +
  scale_colour_manual(values =  RColorBrewer::brewer.pal(9, "Paired")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = -40, hjust = 0)) +
  facet_wrap(vars(metric))


##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
# 
# 
# classification_report %>% 
#   pivot_longer(cols = precision:f1.score, names_to = "metric") %>% 
#   # filter(metric == "precision") %>% 
#   arrange(desc(value)) %>% 
#   ggplot(aes( y = value,x=metric, color = method, label = round(value,3) ) ) +
#   geom_point(alpha = 0.8, size = 2.5,) + 
#   geom_text(size = 2,nudge_x = 0.5) +
#   facet_wrap(vars(X))

















# 
# pheatmap::pheatmap(table(symphony_res$meta_data$label, symphony_res$meta_data$cell_type_pred_knn),
#                    cluster_rows = F,
#                    cluster_cols = F,
#                    scale = 'column',
#                    fontsize = 6)
# 
# 
# 
# 
# pheatmap::pheatmap(table(scmapcell_res$test_labels, scmapcell_res$prediciton$pred_labels_cell),
#                    cluster_rows = F,
#                    cluster_cols = F,
#                    scale = 'column',
#                    fontsize = 6)
# 
# 
# 
# pheatmap::pheatmap(table(scmapcell_res$test_labels, singler_res$prediciton$pred_labels$labels),
#                    cluster_rows = F,
#                    cluster_cols = F,
#                    scale = 'column',
#                    fontsize = 6)
# 
# 
# # Test if order is the same for test labels
# all(scmapcell_res$test_labels == symphony_res$meta_data$label)
# 
# pred_summary <- data.frame("test_label"     = symphony_res$meta_data$label,
#                            "symphony_pred"  = symphony_res$meta_data$cell_type_pred_knn,
#                            "scmapcell_pred" = scmapcell_res$prediciton$pred_labels_cell,
#                            "singler_pred"   = singler_res$prediciton$pred_labels$labels
#                            )
# 
# write.csv(pred_summary, file = "/data/gruen/herman/scDatasets/run_CellAnnotation_tools/prediction_summary.csv")
# 





