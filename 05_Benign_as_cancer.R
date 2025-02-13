#This code will reproduce the analysis when benign is considered cancer/benign class and comparisons
#Fig-6A,6B,6C,6D,6E,6F,6G,6H can be generated.

# BenignAsCancer ----------------------------------------------------------

# GSE106817_Neo -----------------------------------------------------------

meta_GSE106817_neo <- read.csv("neo/GSE106817_metatestfull.csv", stringsAsFactors = TRUE, row.names = 1 )
phenoData_GSE106817_neo <- new("AnnotatedDataFrame",data=meta_GSE106817_neo)

GSE106817_expset_neo <-ExpressionSet(GSE106817mat, phenoData = phenoData_GSE106817_neo)

results_GSE106817_neo <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                                 Data = GSE106817_expset_neo,
                                                 tolerate_missed_genes = TRUE,
                                                 weighted_votes = TRUE,
                                                 verbose = TRUE)
knitr::kable(head(results_GSE106817_neo))
confusion_GSE106817_neo <- caret::confusionMatrix(data = factor(results_GSE106817_neo$max_score, 
                                                                levels = unique(train_object$data$Labels)),
                                               reference = factor(pData(GSE106817_expset_neo)[,"class"],levels = unique(train_object$data$Labels)),mode="everything")
print(confusion_GSE106817_neo)
expr_subset_GSE106817_neo <- GSE106817mat[rownames(GSE106817mat) %in% miRNAs_to_plot, ]

# GSE112264_Neo -----------------------------------------------------------

meta_GSE112264_neo <- read.csv("neo/GSE112264_metatest.csv", stringsAsFactors = TRUE, row.names = 1 )
phenoData_GSE112264_neo <- new("AnnotatedDataFrame",data=meta_GSE112264_neo)

GSE112264_expset_neo <-ExpressionSet(GSE112264mat, phenoData = phenoData_GSE112264_neo)

results_GSE112264_neo <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                                 Data = GSE112264_expset_neo,
                                                 tolerate_missed_genes = TRUE,
                                                 weighted_votes = TRUE,
                                                 verbose = TRUE)
knitr::kable(head(results_GSE112264_neo))

confusion_GSE112264_neo <- caret::confusionMatrix(data = factor(results_GSE112264_neo$max_score, 
                                                                levels = unique(train_object$data$Labels)),
                                                  reference = factor(pData(GSE112264_expset_neo)[,"class"], 
                                                                     levels = unique(train_object$data$Labels)),
                                                  mode="everything")
print(confusion_GSE112264_neo)

plot_confusion_matrix(confusion_GSE112264_neo, "Confusion Matrix for GSE112264_neo")

expr_subset_GSE112264_neo <- GSE112264mat[rownames(GSE112264mat) %in% miRNAs_to_plot, ]

# GSE113486_Neo -----------------------------------------------------------

meta_GSE113486_neo <- read.csv("neo/GSE113486_metatest.csv", stringsAsFactors = TRUE, row.names = 1 )
meta_GSE113486_neo <- meta_GSE113486_h
meta_GSE113486_neo$description <- NULL
phenoData_GSE113486_neo <- new("AnnotatedDataFrame",data=meta_GSE113486_h)

GSE113486_expset_neo <-ExpressionSet(GSE113486mat, phenoData = phenoData_GSE113486_neo)

results_GSE113486_neo <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                                 Data = GSE113486_expset_neo,
                                                 tolerate_missed_genes = TRUE,
                                                 weighted_votes = TRUE,
                                                 verbose = TRUE)
knitr::kable(head(results_GSE113486_neo))

confusion_GSE113486_neo <- caret::confusionMatrix(data = factor(results_GSE113486_neo$max_score, 
                                                                levels = unique(train_object$data$Labels)),
                                                  reference = factor(pData(GSE113486_expset_neo)[,"class"], 
                                                                     levels = unique(train_object$data$Labels)),
                                                  mode="everything")
print(confusion_GSE113486_neo)

plot_confusion_matrix(confusion_GSE113486_neo, "Confusion Matrix for GSE113486_neo")

expr_subset_GSE113486_neo <- GSE113486mat[rownames(GSE113486mat) %in% miRNAs_to_plot, ]

# GSE113740_Neo -----------------------------------------------------------

meta_GSE113740_neo <- read.csv("neo/GSE113740_metatest.csv", stringsAsFactors = TRUE, row.names = 1 )
phenoData_GSE113740_neo <- new("AnnotatedDataFrame",data=meta_GSE113740_neo)

GSE113740_expset_neo <-ExpressionSet(GSE113740mat, phenoData = phenoData_GSE113740_neo)

results_GSE113740_neo <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                                 Data = GSE113740_expset_neo,
                                                 tolerate_missed_genes = TRUE,
                                                 weighted_votes = TRUE,
                                                 verbose = TRUE)
knitr::kable(head(results_GSE113740_neo))

confusion_GSE113740_neo <- caret::confusionMatrix(data = factor(results_GSE113740_neo$max_score, 
                                                                levels = unique(train_object$data$Labels)),
                                                  reference = factor(pData(GSE113740_expset_neo)[,"class"], 
                                                                     levels = unique(train_object$data$Labels)),
                                                  mode="everything")
print(confusion_GSE113740_neo)


plot_confusion_matrix(confusion_GSE113740_neo, "Confusion Matrix for GSE113740_neo")

ggsave("neoplasm/confusion_GSE113740_neo.png")

expr_subset_GSE113740_neo <- GSE113740mat[rownames(GSE113740mat) %in% miRNAs_to_plot, ]


# GSE122497_neo -----------------------------------------------------------

meta_GSE122497_neo <- read.csv("neo/GSE122497_meta.csv", stringsAsFactors = TRUE, row.names = 1 )
phenoData_GSE122497_neo <- new("AnnotatedDataFrame",data=meta_GSE122497_neo)

GSE122497_expset_neo <-ExpressionSet(GSE122497mat, phenoData = phenoData_GSE122497_neo)

results_GSE122497_neo <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                                 Data = GSE122497_expset_neo,
                                                 tolerate_missed_genes = TRUE,
                                                 weighted_votes = TRUE,
                                                 verbose = TRUE)
knitr::kable(head(results_GSE122497_neo))

confusion_GSE122497_neo <- caret::confusionMatrix(data = factor(results_GSE122497_neo$max_score, 
                                                                levels = unique(train_object$data$Labels)),
                                                  reference = factor(pData(GSE122497_expset_neo)[,"class"], 
                                                                     levels = unique(train_object$data$Labels)),
                                                  mode="everything")
print(confusion_GSE122497_neo)

plot_confusion_matrix(confusion_GSE122497_neo, "Confusion Matrix for GSE122497_neo")
expr_subset_GSE122497_neo <- GSE122497mat[rownames(GSE122497mat) %in% miRNAs_to_plot, ]


# GSE137140_Neo -----------------------------------------------------------

meta_GSE137140_neo <- read.csv("neo/GSE137140_metatest.csv", stringsAsFactors = TRUE, row.names = 1 )
phenoData_GSE137140_neo <- new("AnnotatedDataFrame",data=meta_GSE137140_neo)

GSE137140_expset_neo <-ExpressionSet(GSE137140mat, phenoData = phenoData_GSE137140_neo)

results_GSE137140_neo <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                                 Data = GSE137140_expset_neo,
                                                 tolerate_missed_genes = TRUE,
                                                 weighted_votes = TRUE,
                                                 verbose = TRUE)
knitr::kable(head(results_GSE137140_neo))

confusion_GSE137140_neo <- caret::confusionMatrix(data = factor(results_GSE137140_neo$max_score, 
                                                                levels = unique(train_object$data$Labels)),
                                                  reference = factor(pData(GSE137140_expset_neo)[,"class"], 
                                                                     levels = unique(train_object$data$Labels)),
                                                  mode="everything")
print(confusion_GSE137140_neo)


plot_confusion_matrix(confusion_GSE137140_neo, "Confusion Matrix for GSE137140_neo")
ggsave("neoplasm/confusion_GSE137140_neo.png")

expr_subset_GSE137140_neo <- GSE137140mat[rownames(GSE137140mat) %in% miRNAs_to_plot, ]


# GSE139031_Neo -----------------------------------------------------------

meta_GSE139031_neo <- read.csv("neo/GSE139031_metatestall.csv", stringsAsFactors = TRUE, row.names = 1 )
phenoData_GSE139031_neo <- new("AnnotatedDataFrame",data=meta_GSE139031_neo)

GSE139031_expset_neo <-ExpressionSet(GSE139031mat, phenoData = phenoData_GSE139031_neo)

results_GSE139031_neo <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                                 Data = GSE139031_expset_neo,
                                                 tolerate_missed_genes = TRUE,
                                                 weighted_votes = TRUE,
                                                 verbose = TRUE)
knitr::kable(head(results_GSE139031_neo))

confusion_GSE139031_neo <- caret::confusionMatrix(data = factor(results_GSE139031_neo$max_score, 
                                                                levels = unique(train_object$data$Labels)),
                                                  reference = factor(pData(GSE139031_expset_neo)[,"class"], 
                                                                     levels = unique(train_object$data$Labels)),
                                                  mode="everything")
print(confusion_GSE139031_neo)

plot_confusion_matrix(confusion_GSE139031_neo, "Confusion Matrix for GSE139031_neo")
ggsave("neoplasm/confusion_GSE139031_neo.png")

expr_subset_GSE139031_neo <- GSE139031mat[rownames(GSE139031mat) %in% miRNAs_to_plot, ]


# GSE164174_Neo -----------------------------------------------------------

meta_GSE164174_neo <- read.csv("neo/GSE164174_metatest.csv", stringsAsFactors = TRUE, row.names = 1 )
phenoData_GSE164174_neo <- new("AnnotatedDataFrame",data=meta_GSE164174_neo)

GSE164174_expset_neo <-ExpressionSet(GSE164174mat, phenoData = phenoData_GSE164174_neo)

results_GSE164174_neo <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                                 Data = GSE164174_expset_neo,
                                                 tolerate_missed_genes = TRUE,
                                                 weighted_votes = TRUE,
                                                 verbose = TRUE)
knitr::kable(head(results_GSE164174_neo))

confusion_GSE164174_neo <- caret::confusionMatrix(data = factor(results_GSE164174_neo$max_score, 
                                                                levels = unique(train_object$data$Labels)),
                                                  reference = factor(pData(GSE164174_expset_neo)[,"class"], 
                                                                     levels = unique(train_object$data$Labels)),
                                                  mode="everything")
print(confusion_GSE164174_neo)

plot_confusion_matrix(confusion_GSE164174_neo, "Confusion Matrix for GSE164174_neo")
expr_subset_GSE164174_neo <- GSE164174mat[rownames(GSE164174mat) %in% miRNAs_to_plot, ]

# GSE73002_Neo ------------------------------------------------------------

meta_GSE73002_neo <- read.csv("neo/metadata73002_dup_rem.csv", stringsAsFactors = TRUE, row.names = 1 )
phenoData_GSE73002_neo <- new("AnnotatedDataFrame",data=meta_GSE73002_neo)

GSE73002_expset_neo <-ExpressionSet(GSE73002mat, phenoData = phenoData_GSE73002_neo)

results_GSE73002_neo <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                                Data = GSE73002_expset_neo,
                                                tolerate_missed_genes = TRUE,
                                                weighted_votes = TRUE,
                                                verbose = TRUE)
knitr::kable(head(results_GSE73002_neo))

confusion_GSE73002_neo <- caret::confusionMatrix(data = factor(results_GSE73002_neo$max_score, 
                                                               levels = unique(train_object$data$Labels)),
                                                 reference = factor(pData(GSE73002_expset_neo)[,"class"], 
                                                                    levels = unique(train_object$data$Labels)),
                                                 mode="everything")
print(confusion_GSE73002_neo)

plot_confusion_matrix(confusion_GSE73002_neo, "Confusion Matrix for GSE73002_neo")
expr_subset_GSE73002_neo <- GSE73002mat[rownames(GSE73002mat) %in% miRNAs_to_plot, ]




# 6A_Benignasnon_cancer_test ----------------------------------------------

ggplot(all_cm, aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), color = "black", fontface = "bold", size = 14) +  # Bold the text in tiles
  facet_wrap(~ Dataset, ncol = 3) +
  scale_fill_gradient(low = "lavender", high = "lightblue") +
  theme_minimal() +
  labs(title = "Confusion Matrices - Test Datasets",
       x = "Predicted Label", y = "True Label") +
  theme(
    plot.title = element_text(face = "bold", size = 6, hjust = 0.5),  # Bold and center the title
    axis.title.x = element_text(face = "bold", size = 24, color = "black"),  # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 24, color = "black"),  # Bold y-axis title
    axis.text.x = element_text(face = "bold", size = 22, color = "black"),   # Bold x-axis text
    axis.text.y = element_text(face = "bold", size = 22, color = "black"),   # Bold y-axis text
    strip.text = element_text(face = "bold", size = 28, color = "black"),    # Bold facet titles
    legend.text = element_text(face = "bold", size = 18),  # Bold and increase legend text size
    legend.title = element_text(face = "bold", size = 16)  # Bold and increase legend title size
  )

ggsave("Figures/Confusion_benignasnon-cancer_test.png")

# 6B_ConfusionMatrix_BenignasCancer ----------------------------------------------

confusion_matrices_neo <- list(
  GSE122497 = confusion_GSE122497_neo,
  GSE106817 = confusion_GSE106817_neo,
  GSE137140 = confusion_GSE137140_neo,
  GSE164174 = confusion_GSE164174_neo,
  GSE113740 = confusion_GSE113740_neo,
  GSE112264 = confusion_GSE112264_neo,
  GSE113486 = confusion_GSE113486_neo,
  GSE139031 = confusion_GSE139031_neo,
  GSE73002 = confusion_GSE73002_neo)

# Convert each matrix to a data frame and add a `Dataset` column
all_cm_neo <- do.call(rbind, lapply(names(confusion_matrices_neo), function(name) {
  cm_neo <- as.data.frame(confusion_matrices_neo[[name]]$table)  # assumes confusion_matrix$table format
  cm_neo$Dataset <- name  # add a column to indicate the dataset
  cm_neo  # return modified data frame
}))

all_cm_neo <- all_cm_neo %>%
  mutate(
    Prediction = recode(Prediction, "cancer" = "cancer/benign", "non_cancer" = "non-cancer/non-benign"),
    Reference = recode(Reference, "cancer" = "cancer/benign", "non_cancer" = "non-cancer/non-benign")
  )

# Now plot with ggplot
ggplot(all_cm_neo, aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), color = "black", fontface = "bold", size = 22) +  # Bold the text in tiles
  facet_wrap(~ Dataset, ncol = 3) +
  scale_fill_gradient(low = "lavender", high = "lightblue") +
  theme_minimal() +
  labs(title = "Confusion Matrices - Test Datasets with Benign as cancer",
       x = "Predicted Label", y = "True Label") +
  theme(
    plot.title = element_text(face = "bold", size = 4, hjust = 0.5),  # Bold and center the title
    axis.title.x = element_text(face = "bold", size = 32, color = "black"),  # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 32, color = "black"),  # Bold y-axis title
    axis.text.x = element_text(face = "bold", size = 18, color = "black"),   # Bold x-axis text
    axis.text.y = element_text(face = "bold", size = 30, color = "black"),   # Bold y-axis text
    strip.text = element_text(face = "bold", size = 36, color = "black"),    # Bold facet titles
    legend.text = element_text(face = "bold", size = 18),  # Bold and increase legend text size
    legend.title = element_text(face = "bold", size = 26)  # Bold and increase legend title size
  )

ggsave("Figures/Confusion_Matrix_withbenignascancer_testlabelsdata.png", dpi = 600, bg = "white", height = 18, width = 24)



# 6C_Train_benignas_non_cancer --------------------------------------------

plot_confusion_matrix(confusion_train, "Confusion Matrix for Training data")
ggsave("Figures/Confusion_benignasnon-cancer_train.png")


# 6D_Train_benign_Visualization ----------------------------------------------

metadata_train_benign <- read.csv("neo/metadata_full_benign_marked.csv", row.names = 1, stringsAsFactors = TRUE)
umap_coords_benign <- as.data.frame(umap_result$layout)  # Extract coordinates
colnames(umap_coords_benign) <- c("UMAP1", "UMAP2")
umap_coords_benign$State <- metadata_train_benign$class[match(rownames(umap_coords_benign), metadata.subset_GSE211692$geo_accession)]
umap_coords_benign$Category <- ifelse(umap_coords_benign$State == "cancer", "cancer",
                                      ifelse(umap_coords_benign$State == "benign", "benign", "non_cancer"))

tsne_coords_benign <- as.data.frame(tsne_result$Y)
tsne_coords_benign$Category <- umap_coords_benign$Category  # Add cancer/non-cancer/benign labels

# Plot t-SNE
ggplot(tsne_coords_benign, aes(x = V1, y = V2, color = Category)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Plot of Cancer, Non-Cancer, and Benign",
       x = "t-SNE 1",
       y = "t-SNE 2") +
  scale_color_manual(values = c("cancer" = "red", "non_cancer" = "blue", "benign" = "darkorange")) +
  theme(legend.position = "right")

ggsave("Figures/t-sne_train_with_benign.png")


#Train_data_when benign is considered "cancer/benign" group

meta_neo <- read.csv("neo/meta_train_full_flip.csv", row.names = 1, stringsAsFactors = TRUE)
phenoData_neo<- new("AnnotatedDataFrame",data=meta_neo)

neoset <- ExpressionSet(highly_expressed_train, phenoData = phenoData_neo)

results_neotrain <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                            Data = neoset,
                                            tolerate_missed_genes = TRUE,
                                            weighted_votes = TRUE,
                                            verbose = TRUE)
knitr::kable(head(results_neotrain))

confusion_neotrain <- caret::confusionMatrix(data = factor(results_neotrain$max_score, 
                                                           levels = unique(train_object$data$Labels)),
                                             reference = factor(pData(neoset)[,"class"], 
                                                                levels = unique(train_object$data$Labels)),
                                             mode="everything")
print(confusion_neotrain)

cm_table_neotrain <- as.data.frame(confusion_neotrain$table)
cm_table_neotrain$Prediction <- factor(cm_table_neotrain$Prediction, 
                                       levels = unique(cm_table_neotrain$Prediction))
levels(cm_table_neotrain$Prediction) <- c("cancer/benign", "non_cancer/non_benign")

cm_table_neotrain$Reference <- factor(cm_table_neotrain$Reference, 
                                      levels = unique(cm_table_neotrain$Reference))
levels(cm_table_neotrain$Reference) <- c("cancer/benign", "non_cancer/non_benign")

ggplot(cm_table_neotrain, aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), color = "black", fontface = "bold", size = 18) + # Display frequency in tiles
  scale_fill_gradient(low = "lavender", high = "lightblue") +
  theme_minimal() +
  labs(title = "Confusion Matrix for Neo Datasets",
       x = "Predicted Label", y = "True Label") +
  theme(
    plot.title = element_text(face = "bold", size = 8, hjust = 0.5),  # Bold and center the title
    axis.title.x = element_text(face = "bold", size = 32, color = "black"),  # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 32, color = "black"),  # Bold y-axis title
    axis.text.x = element_text(face = "bold", size = 28, color = "black"),   # Bold x-axis text
    axis.text.y = element_text(face = "bold", size = 24, color = "black"),   # Bold y-axis text
    strip.text = element_text(face = "bold", size = 20, color = "black"),    # Bold facet titles
    legend.text = element_text(face = "bold", size = 16),  # Bold and increase legend text size
    legend.title = element_text(face = "bold", size = 16)  # Bold and increase legend title size
  )
ggsave("merged_confusion_Train_benign as cancer.png")


# 6E_Test_Benignasnon_cancer_merged_conf ----------------------------------

plot_confusion_matrix(confusion_test_merged, "Confusion Matrix for Test_Datasets")
ggsave("Figures/confusion_merged_Test.png")



# 6F_ Merged_confusion_benignascancer_merged_conf --------------------------------------------------------

# List of neo datasets
datasets_neo <- list(
  GSE106817 = GSE106817_expset_neo,
  GSE112264 = GSE112264_expset_neo,
  GSE113486 = GSE113486_expset_neo,
  GSE113740 = GSE113740_expset_neo,
  GSE122497 = GSE122497_expset_neo,
  GSE137140 = GSE137140_expset_neo,
  GSE139031 = GSE139031_expset_neo,
  GSE164174 = GSE164174_expset_neo,
  GSE73002 = GSE73002_expset_neo
)

common_features_neo <- Reduce(intersect, lapply(datasets_neo, function(x) featureNames(x)))
datasets_neo_subset <- lapply(datasets_neo, function(x) x[common_features_neo, ])
expression_matrices_neo <- lapply(datasets_neo_subset, exprs)
merged_exprs_neo <- do.call(cbind, expression_matrices_neo)

metadata_list_neo <- lapply(names(datasets_neo_subset), function(dataset_name) {
  meta <- pData(datasets_neo_subset[[dataset_name]])
  meta$Dataset <- dataset_name  # Add dataset name as a new column
  return(meta)
})

metadata_list_neo[[3]] <- metadata_list_neo[[3]][, !(colnames(metadata_list_neo[[3]]) %in% c("description"))]
merged_metadata_neo <- do.call(rbind, metadata_list_neo)
merged_test_mat_neo <- as.matrix(merged_exprs_neo)

phenoData_mergedtest_neo <- new("AnnotatedDataFrame", data = merged_metadata_neo)
merged_expressionset_neo <- ExpressionSet(merged_test_mat_neo, phenoData = phenoData_mergedtest_neo)

#Predict using the classifier on the merged neo expression set
results_test_neo <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                            Data = merged_expressionset_neo,
                                            tolerate_missed_genes = TRUE,
                                            weighted_votes = TRUE,
                                            verbose = TRUE)
knitr::kable(head(results_test_neo))

#ConfusionMatrix

confusion_test_neo <- caret::confusionMatrix(
  data = factor(results_test_neo$max_score, 
                levels = unique(train_object$data$Labels)),
  reference = factor(pData(merged_expressionset_neo)[, "class"], 
                     levels = unique(train_object$data$Labels)),
  mode = "everything"
)
print(confusion_test_neo)
plot_confusion_matrix(confusion_test_neo, "Confusion Matrix for Neo Datasets")

ggsave("merged_confusion_Test_benign as cancer.png")

# 6G_PerformanceMetrices_BenignasCancer_Test --------------------------------------

cm_metrics_neo <- data.frame(
  Dataset = names(confusion_matrices_neo),
  Accuracy = sapply(confusion_matrices_neo, function(cm) {
    accuracy <- sum(diag(cm$table)) / sum(cm$table)
    return(accuracy)
  }),
  Precision = sapply(confusion_matrices_neo, function(cm) {
    precision <- cm$table[2, 2] / sum(cm$table[2, ])
    return(precision)
  }),
  Recall = sapply(confusion_matrices_neo, function(cm) {
    recall <- cm$table[2, 2] / sum(cm$table[, 2])
    return(recall)
  }),
  F1 = sapply(confusion_matrices_neo, function(cm) {
    precision <- cm$table[2, 2] / sum(cm$table[2, ])
    recall <- cm$table[2, 2] / sum(cm$table[, 2])
    f1 <- 2 * (precision * recall) / (precision + recall)
    return(f1)
  }),
  Sensitivity = sapply(confusion_matrices_neo, function(cm) {
    sensitivity <- cm$table[2, 2] / sum(cm$table[, 2])  # True Positive Rate
    return(sensitivity)
  }),
  Specificity = sapply(confusion_matrices_neo, function(cm) {
    specificity <- cm$table[1, 1] / sum(cm$table[1, ])  # True Negative Rate
    return(specificity)
  })
)

cm_metrics_long_neo <- pivot_longer(cm_metrics_neo, cols = -Dataset, names_to = "Metric", values_to = "Value")
cm_metrics_long_neo$Dataset <- factor(cm_metrics_long_neo$Dataset, 
                                      levels = rev(c("GSE106817", "GSE112264", "GSE113486", 
                                                     "GSE113740", "GSE122497", "GSE137140", 
                                                     "GSE139031", "GSE164174", "GSE73002")))

print(levels(cm_metrics_long_neo$Dataset))

ggplot(cm_metrics_long_neo, aes(x = Metric, y = Dataset, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = percent(Value, accuracy = 0.1)), 
            color = "black", size = 5, fontface = "bold") +  # Bold percentage labels
  scale_fill_gradient(low = "lavender", high = "light blue") +  # Custom gradient from light blue to darker blue
  labs(title = "Rules-based Classifier Metrics - Test Datasets with Benign as Cancer", 
       x = "Metric", 
       y = "Dataset", 
       fill = "Metric Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14, color = "black"),  # Bold x-axis labels
    axis.text.y = element_text(face = "bold", size = 14, color = "black"),                         # Bold y-axis labels
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),             # Bold and centered title
    axis.title.x = element_text(size = 18, face = "bold", color = "black"),                        # Bold x-axis title
    axis.title.y = element_text(size = 18, face = "bold", color = "black"),                        # Bold y-axis title
    strip.text = element_text(size = 12, face = "bold"),                          # Bold facet labels
    legend.text = element_text(size = 10, face = "bold"),                         # Bold legend text
    legend.title = element_text(size = 12, face = "bold")                         # Bold legend title
  ) +
  scale_x_discrete(expand = c(0.1, 0.1)) +  # Add some spacing for x-axis labels
  scale_y_discrete(limits = sort(unique(cm_metrics_long$Dataset)))  # Sort y-axis labels in ascending order

ggsave("Figures/wb_heatmap_metrices_benignneo.png")


# 6H_ROC_BenignasCancer ------------------------------------------------------

true_labels_122497_neo <- as.numeric(factor(meta_GSE122497_neo$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_122497_neo <- results_GSE122497_neo$cancer
roc_GSE122497_neo <- roc(true_labels_122497_neo, predicted_probs_122497_neo)
plot(roc_GSE122497_neo, col = "blue", main = "ROC Curve for GSE122497", lwd = 2)
auc_value_122497_neo <- auc(roc_GSE122497_neo)
text(0.6, 0.2, paste("AUC =", round(auc_value_122497_neo, 4)), col = "blue", cex = 1.2)

true_labels_106817_neo <- as.numeric(factor(meta_GSE106817_neo$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_106817_neo <- results_GSE106817_neo$cancer
roc_GSE106817_neo <- roc(true_labels_106817_neo, predicted_probs_106817_neo)
plot(roc_GSE106817_neo, col = "blue", main = "ROC Curve for GSE106817", lwd = 2)
auc_value_106817_neo <- auc(roc_GSE106817_neo)
text(0.6, 0.2, paste("AUC =", round(auc_value_106817_neo, 4)), col = "blue", cex = 1.2)


true_labels_137140_neo <- as.numeric(factor(meta_GSE137140_neo$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_137140_neo <- results_GSE137140_neo$cancer
roc_GSE137140_neo <- roc(true_labels_137140_neo, predicted_probs_137140_neo)
plot(roc_GSE137140_neo, col = "blue", main = "ROC Curve for GSE137140", lwd = 2)
auc_value_137140_neo <- auc(roc_GSE137140_neo)
text(0.6, 0.2, paste("AUC =", round(auc_value_137140_neo, 4)), col = "blue", cex = 1.2)

true_labels_164174_neo <- as.numeric(factor(meta_GSE164174_neo$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_164174_neo <- results_GSE164174_neo$cancer
roc_GSE164174_neo <- roc(true_labels_164174_neo, predicted_probs_164174_neo)
plot(roc_GSE164174_neo, col = "blue", main = "ROC Curve for GSE164174", lwd = 2)
auc_value_164174_neo <- auc(roc_GSE164174_neo)
text(0.6, 0.2, paste("AUC =", round(auc_value_164174_neo, 4)), col = "blue", cex = 1.2)

true_labels_113740_neo <- as.numeric(factor(meta_GSE113740_neo$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_113740_neo <- results_GSE113740_neo$cancer
roc_GSE113740_neo <- roc(true_labels_113740_neo, predicted_probs_113740_neo)
plot(roc_GSE113740_neo, col = "blue", main = "ROC Curve for GSE113740", lwd = 2)
auc_value_113740_neo <- auc(roc_GSE113740_neo)
text(0.6, 0.2, paste("AUC =", round(auc_value_113740_neo, 4)), col = "blue", cex = 1.2)


true_labels_112264_neo <- as.numeric(factor(meta_GSE112264_neo$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_112264_neo <- results_GSE112264_neo$cancer
roc_GSE112264_neo <- roc(true_labels_112264_neo, predicted_probs_112264_neo)
plot(roc_GSE112264_neo, col = "blue", main = "ROC Curve for GSE112264", lwd = 2)
auc_value_112264_neo <- auc(roc_GSE112264_neo)
text(0.6, 0.2, paste("AUC =", round(auc_value_112264_neo, 4)), col = "blue", cex = 1.2)


true_labels_113486_neo <- as.numeric(factor(meta_GSE113486_neo$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_113486_neo <- results_GSE113486_neo$cancer
roc_GSE113486_neo <- roc(true_labels_113486_neo, predicted_probs_113486_neo)
plot(roc_GSE113486_neo, col = "blue", main = "ROC Curve for GSE113486", lwd = 2)
auc_value_113486_neo <- auc(roc_GSE113486_neo)
text(0.6, 0.2, paste("AUC =", round(auc_value_113486_neo, 4)), col = "blue", cex = 1.2)

true_labels_139031_neo <- as.numeric(factor(meta_GSE139031_neo$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_139031_neo <- results_GSE139031_neo$cancer
roc_GSE139031_neo <- roc(true_labels_139031_neo, predicted_probs_139031_neo)
plot(roc_GSE139031_neo, col = "blue", main = "ROC Curve for GSE139031", lwd = 2)
auc_value_139031_neo <- auc(roc_GSE139031)
text(0.6, 0.2, paste("AUC =", round(auc_value_139031_neo, 4)), col = "blue", cex = 1.2)


true_labels_73002_neo <- as.numeric(factor(meta_GSE73002_neo$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_73002_neo <- results_GSE73002_neo$cancer
roc_GSE73002_neo <- roc(true_labels_73002_neo, predicted_probs_73002_neo)
plot(roc_GSE73002_neo, col = "blue", main = "ROC Curve for GSE73002", lwd = 2)
auc_value_73002_neo <- auc(roc_GSE73002_neo)
text(0.6, 0.2, paste("AUC =", round(auc_value_73002_neo, 4)), col = "blue", cex = 1.2)


roc_data_neo <- list(
  GSE122497 = roc_GSE122497_neo,
  GSE106817 = roc_GSE106817_neo,
  GSE137140 = roc_GSE137140_neo,
  GSE164174 = roc_GSE164174_neo,
  GSE113740 = roc_GSE113740_neo,
  GSE112264 = roc_GSE112264_neo,
  GSE113486 = roc_GSE113486_neo,
  GSE139031 = roc_GSE139031_neo,
  GSE73002 = roc_GSE73002_neo
)

# Convert each ROC object to a data frame and add a `Dataset` column
all_roc_neo <- do.call(rbind, lapply(names(roc_data_neo), function(name) {
  roc_neo <- roc_data_neo[[name]]
  df <- data.frame(
    FPR = 1 - roc_neo$specificities,  # False Positive Rate
    TPR = roc_neo$sensitivities,      # True Positive Rate
    Threshold = roc_neo$thresholds    # Thresholds for ROC
  )
  df$Dataset <- name  # add dataset name as a column
  
  # Calculate AUC for each ROC object
  auc_value_neo <- auc(roc_neo)  # Compute the AUC for the current ROC curve
  
  df$AUC <- auc_value_neo  # Add the AUC value to the dataframe for later use
  return(df)
}))


ggplot(all_roc_neo, aes(x = FPR, y = TPR, color = Dataset)) +
  geom_line() +
  
  # Add AUC value near the top-right corner (adjusted for visibility)
  geom_label(data = subset(all_roc_neo, FPR == max(FPR) & TPR == max(TPR)), 
             aes(label = paste("AUC = ", format(AUC, digits = 4))),
             color = "white", size = 2.5, fontface = "bold", 
             fill = "black", label.padding = unit(0.5, "lines"),
             label.r = unit(0.15, "lines"), hjust = 1.5, vjust = 1.5) +
  
  # Add Dataset name in large font at the top-left corner
  geom_text(data = subset(all_roc_neo, FPR == min(FPR) & TPR == max(TPR)),
            aes(label = Dataset), color = "black", size = 7, fontface = "bold", 
            hjust = -0.2, vjust = 1.5) + 
  
  theme_minimal() +
  labs(title = "ROC Curves for Multiple Datasets",
       x = "False Positive Rate", y = "True Positive Rate") +
  scale_color_viridis_d() +
  facet_wrap(~ Dataset, ncol = 3) +  # Create a 3-column grid layout
  theme(legend.position = "none",  # Remove legend as we're using color for dataset
        strip.text = element_text(size = 12),  # Adjust text size for dataset labels
        plot.margin = margin(5, 5, 5, 5))  # Adjust plot margins



# Precompute AUC values and positions for labels
auc_positions_neo <- all_roc_neo %>%
  group_by(Dataset) %>%
  summarise(
    FPR_AUC = max(FPR), TPR_AUC = max(TPR),  # Position for AUC label
    AUC_Value = unique(AUC)  # Extract the AUC value for the dataset
  ) %>%
  arrange(Dataset)  # Ensure the datasets are in order

all_roc_neo$Dataset <- factor(all_roc_neo$Dataset, levels = auc_positions_neo$Dataset)

# Plot
ggplot(all_roc_neo, aes(x = FPR, y = TPR, color = Dataset)) +
  geom_line(size = 1) +
  
  # Add AUC value near the top-right corner for each dataset
  geom_label(data = auc_positions_neo,
             aes(x = FPR_AUC, y = TPR_AUC,
                 label = paste("AUC = ", format(AUC_Value, digits = 4))),
             inherit.aes = FALSE,  # Prevent aesthetics inheritance
             color = "white", fill = "black", size = 5, fontface = "bold",
             label.padding = unit(0.5, "lines"), hjust = 1, vjust = 1) +
  
  theme_minimal() +
  labs(title = "ROC Curves for Test Datasets with benign as cancer",
       x = "False Positive Rate", y = "True Positive Rate") +
  scale_color_viridis_d() +
  facet_wrap(~ Dataset, ncol = 3) +  # Create a 3-column grid layout
  theme(
    legend.position = "none",  # Remove legend as we're using color for dataset
    strip.text = element_text(size = 16, face = "bold"),  # Bold facet labels
    axis.title.x = element_text(size = 18, face = "bold"),  # Bold x-axis title
    axis.title.y = element_text(size = 18, face = "bold"),  # Bold y-axis title
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),  # Bold x-axis text
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),  # Bold y-axis text
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Bold and centered title
    plot.margin = margin(5, 5, 5, 5)  # Adjust plot margins
  )

ggsave("roc_plots_bold_benignascancer.png", dpi = 600, bg = "white", height = 8, width = 16)

