#This code will reproduce the pre-processing of Test Data and the classifier performance
#Fig-4A,4B,4C,4D,4E can be generated.
########################Preprocessing###########################################

# Test_datasets -----------------------------------------------------------

# GSE106817_Test ----------------------------------------------------------

# Fetch the GEO data
gse106817 <- getGEO("GSE106817", GSEMatrix = TRUE)
metadata_106817 = pData(phenoData(gse106817[[1]]))
metadata_106817.subset = dplyr::select(metadata_106817,c(1,19))

# organize
metadata.subset_cancer = metadata_106817 %>%
  dplyr::select(1,19) %>%
  dplyr::rename(state = description) %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))
str(gse106817)

gse106817_data <- gse106817[[1]]
meta_GSE106817 <- pData(gse106817_data)
annot_GSE106817 <- fData(gse106817_data)

exp_GSE106817 <- exprs(gse106817_data)
exp_106817 <- as.data.frame(exp_GSE106817)
frequency_table_106817 <- table(metadata_106817.subset$description)
frequency_table_filtered <- table(metadata_106817.subset$description)
frequency_table_filtered
GSE106817_frame <- as.data.frame(exp_106817)
rownames(GSE106817_frame) <- annot_GSE106817$miRNA_ID_LIST
GSE106817_frame <- na.omit(GSE106817_frame)


cleaned_GSE106817 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE106817_frame))
cleaned_GSE106817 <- gsub("-", "_", cleaned_GSE106817) 
cleaned_GSE106817 <- gsub(",", ".", cleaned_GSE106817) 

rownames(GSE106817_frame) <- cleaned_GSE106817
GSE106817mat <- as.matrix(GSE106817_frame)

meta_GSE106817_h <- read.csv("test/GSE106817_metadatah.csv", stringsAsFactors = TRUE, row.names = 1 )
phenoData_GSE106817 <- new("AnnotatedDataFrame",data=meta_GSE106817_h)

#t-SNE Visulaization

umap_result_GSE106817 <- umap(t(GSE106817_frame))# Plot UMAP with specific colors for cancer and non-cancer
umap_coords_GSE106817 <- as.data.frame(umap_result_GSE106817$layout)
colnames(umap_coords_GSE106817) <- c("UMAP1", "UMAP2")
umap_coords_GSE106817$State <- meta_GSE106817_h$class[match(rownames(umap_coords_GSE106817), rownames(meta_GSE106817_h))]
umap_coords_GSE106817$Category <- ifelse(umap_coords_GSE106817$State == "cancer", "cancer", "non_cancer")
unique(umap_coords_GSE106817$State)

tsne_result_GSE106817 <- Rtsne(t(GSE106817_frame), dims = 2, pca = TRUE, check_duplicates = FALSE)
tsne_coords_GSE106817 <- as.data.frame(tsne_result_GSE106817$Y)
tsne_coords_GSE106817$Category <- umap_coords_GSE106817$Category  # Add cancer/non-cancer labels

ggplot(tsne_coords_GSE106817, aes(x = V1, y = V2, color = Category)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Plot GSE106817",
       x = "t-SNE 1",
       y = "t-SNE 2") +
  scale_color_manual(values = c("cancer" = "red", "non_cancer" = "blue")) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12, face = "bold"), # Increase legend text size
    legend.title = element_text(size = 14, face = "bold") # Increase legend title size
  )
ggsave("Figures/t-sne_GSE106817.png")

GSE106817_expset <-ExpressionSet(GSE106817mat, phenoData = phenoData_GSE106817)
results_GSE106817 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                             Data = GSE106817_expset,
                                             tolerate_missed_genes = TRUE,
                                             weighted_votes = TRUE,
                                             verbose = TRUE)
knitr::kable(head(results_GSE106817))

confusion_GSE106817 <- caret::confusionMatrix(data = factor(results_GSE106817$max_score, 
                                                            levels = unique(train_object$data$Labels)),
                                              reference = factor(pData(GSE106817_expset)[,"class"], 
                                                                 levels = unique(train_object$data$Labels)),
                                              mode="everything")
print(confusion_GSE106817)
plot_confusion_matrix(confusion_GSE106817, "Confusion Matrix for GSE106817")

expr_subset_GSE106817 <- GSE106817_frame[rownames(GSE106817_frame) %in% miRNAs_to_plot, ]
expr_long_GSE106817 <- expr_subset_GSE106817 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)

metaplot_GSE106817 <- as.data.frame(meta_GSE106817_h)
metaplot_GSE106817$Sample <- rownames(metaplot_GSE106817)
rownames(metaplot_GSE106817) <- NULL

# Prepare predictions from results
predictions_GSE106817 <- data.frame(
  Sample = rownames(results_GSE106817),  # Use row names as Sample names
  Prediction = results_GSE106817$max_score  # Rename this to match your earlier context
)
head(predictions_GSE106817)

expr_long_GSE106817 <- expr_long_GSE106817 %>%
  left_join(metaplot_GSE106817, by = "Sample") %>%
  left_join(predictions_GSE106817, by = "Sample") %>%
  rename(actual_label = class)

metaplot_GSE106817 <- metaplot_GSE106817 %>%
  rename(actual_label = class)

# Merge the actual labels with predictions
predictions_GSE106817 <- predictions_GSE106817 %>%
  left_join(metaplot_GSE106817, by = "Sample")

prediction_outcomes_GSE106817 <- predictions_GSE106817 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))
head(prediction_outcomes_GSE106817)
table(prediction_outcomes_GSE106817$Correct)


# GSE112264_Test ----------------------------------------------------------

# Fetch the GEO data
gse112264 <- getGEO("GSE112264", GSEMatrix = TRUE)
metadata_112264 = pData(phenoData(gse112264[[1]]))
metadata_112264.subset = dplyr::select(metadata_112264,c(1,40))
frequency_table_112264 <- table(metadata_112264.subset$`disease state:ch1`)

# organize
metadata.subset_cancer = metadata_112264 %>%
  dplyr::select(1,40) %>%
  dplyr::rename(state = "disease state:ch1") %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))
str(gse112264)
gse112264_data <- gse112264[[1]]

# Print the sample information
meta_GSE112264 <- pData(gse112264_data)
annot_GSE112264 <- fData(gse112264_data)
exp_GSE112264 <- exprs(gse112264_data)

GSE112264_frame <- as.data.frame(exp_GSE112264)
rownames(GSE112264_frame) <- annot_GSE112264$miRNA_ID_LIST
GSE112264_frame <- na.omit(GSE112264_frame)

cleaned_GSE112264 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE112264_frame))
cleaned_GSE112264 <- gsub("-", "_", cleaned_GSE112264) 
cleaned_GSE112264 <- gsub(",", ".", cleaned_GSE112264) 

rownames(GSE112264_frame) <- cleaned_GSE112264
GSE112264mat <- as.matrix(GSE112264_frame)

meta_GSE112264_h <- read.csv("test/GSE112264_metadatah.csv", stringsAsFactors = TRUE,row.names = 1)
phenoData_GSE112264 <- new("AnnotatedDataFrame",data=meta_GSE112264_h)

umap_result_GSE112264 <- umap(t(GSE112264_frame))
umap_coords_GSE112264 <- as.data.frame(umap_result_GSE112264$layout)
colnames(umap_coords_GSE112264) <- c("UMAP1", "UMAP2")
umap_coords_GSE112264$State <- meta_GSE112264_h$class[match(rownames(umap_coords_GSE112264), rownames(meta_GSE112264_h))]
umap_coords_GSE112264$Category <- ifelse(umap_coords_GSE112264$State == "cancer", "cancer", "non_cancer")
unique(umap_coords_GSE112264$State)

tsne_result_GSE112264 <- Rtsne(t(GSE112264_frame), dims = 2, pca = TRUE, check_duplicates = FALSE)
tsne_coords_GSE112264 <- as.data.frame(tsne_result_GSE112264$Y)
tsne_coords_GSE112264$Category <- umap_coords_GSE112264$Category  # Add cancer/non-cancer labels
ggplot(tsne_coords_GSE112264, aes(x = V1, y = V2, color = Category)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Plot GSE112264",
       x = "t-SNE 1",
       y = "t-SNE 2") +
  scale_color_manual(values = c("cancer" = "red", "non_cancer" = "blue")) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12, face = "bold"), # Increase legend text size
    legend.title = element_text(size = 14, face = "bold") # Increase legend title size
  )

ggsave("Figures/t-sne_GSE112264.png")

GSE112264_expset <-ExpressionSet(GSE112264mat, phenoData = phenoData_GSE112264)
results_GSE112264 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                             Data = GSE112264_expset,
                                             tolerate_missed_genes = TRUE,
                                             weighted_votes = TRUE,
                                             verbose = TRUE)
knitr::kable(head(results_GSE112264))

confusion_GSE112264 <- caret::confusionMatrix(data = factor(results_GSE112264$max_score, 
                                                            levels = unique(train_object$data$Labels)),
                                              reference = factor(pData(GSE112264_expset)[,"class"], 
                                                                 levels = unique(train_object$data$Labels)),
                                              mode="everything")
print(confusion_GSE112264)
plot_confusion_matrix(confusion_GSE112264, "Confusion Matrix for GSE112264")

expr_subset_GSE112264 <- GSE112264_frame[rownames(GSE112264_frame) %in% miRNAs_to_plot, ]

library(tidyverse)
expr_long_GSE112264 <- expr_subset_GSE112264 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)

metaplot_GSE112264 <- as.data.frame(meta_GSE112264_h)
metaplot_GSE112264$Sample <- rownames(metaplot_GSE112264)
rownames(metaplot_GSE112264) <- NULL  # Optionally, reset the row names

# Prepare predictions from results
predictions_GSE112264 <- data.frame(
  Sample = rownames(results_GSE112264),
  Prediction = results_GSE112264$max_score
)
head(predictions_GSE112264)


expr_long_GSE112264 <- expr_long_GSE112264 %>%
  left_join(metaplot_GSE112264, by = "Sample") %>%
  left_join(predictions_GSE112264, by = "Sample") %>%
  rename(actual_label = class)

metaplot_GSE112264 <- metaplot_GSE112264 %>%
  rename(actual_label = class)

predictions_GSE112264 <- predictions_GSE112264 %>%
  left_join(metaplot_GSE112264, by = "Sample")

prediction_outcomes_GSE112264 <- predictions_GSE112264 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))

head(prediction_outcomes_GSE112264)
table(prediction_outcomes_GSE112264$Correct)

# GSE113486_Test ----------------------------------------------------------

# Fetch the GEO data
GSE113486 <- getGEO("GSE113486", GSEMatrix = TRUE)
metadata_GSE113486 = pData(phenoData(GSE113486[[1]]))
metadata_GSE113486.subset = dplyr::select(metadata_GSE113486,c(1,36))
colnames(metadata_GSE113486.subset) <- c("Title", "Description")
frequency_table_GSE113486 <- table(metadata_GSE113486.subset$Description)

# organize
metadata.subset_cancer = metadata_GSE113486.subset %>%
  dplyr::rename(state = Description) %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))

str(GSE113486)
GSE113486_data <- GSE113486[[1]]

meta_GSE113486 <- pData(GSE113486_data)
annot_GSE113486 <- fData(GSE113486_data)

exp_GSE113486 <- exprs(GSE113486_data)
GSE113486_frame <- as.data.frame(exp_GSE113486)
rownames(GSE113486_frame) <- annot_GSE113486$miRNA_ID_LIST
GSE113486_frame <- na.omit(GSE113486_frame)

#There are duplicate samples identified in this dataset

GSE113486_WO <- read.csv("test/GSE113486_noduplicate.csv")
GSE113486_ID <- GSE113486_WO$ID
GSE11386_framed <- dplyr::select(GSE113486_frame, dplyr::all_of(GSE113486_ID))

cleaned_GSE113486 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE11386_framed))
cleaned_GSE113486 <- gsub("-", "_", cleaned_GSE113486) 
cleaned_GSE113486 <- gsub(",", ".", cleaned_GSE113486) 
rownames(GSE11386_framed) <- cleaned_GSE113486
GSE113486mat <- as.matrix(GSE11386_framed)

metaf_GSE113486 <- read.csv("test/GSE113486_metadatah.csv", row.names = 1, stringsAsFactors = TRUE)
meta_GSE113486_h <- metaf_GSE113486[rownames(metaf_GSE113486) %in% colnames(GSE113486mat), , drop = FALSE]
phenoData_GSE113486 <- new("AnnotatedDataFrame",data=meta_GSE113486_h)

umap_result_GSE113486 <- umap(t(GSE113486mat))
umap_coords_GSE113486 <- as.data.frame(umap_result_GSE113486$layout)
colnames(umap_coords_GSE113486) <- c("UMAP1", "UMAP2")
umap_coords_GSE113486$State <- meta_GSE113486_h$class[match(rownames(umap_coords_GSE113486), rownames(meta_GSE113486_h))]
umap_coords_GSE113486$Category <- ifelse(umap_coords_GSE113486$State == "cancer", "cancer", "non_cancer")
unique(umap_coords_GSE113486$State)

tsne_result_GSE113486 <- Rtsne(t(GSE113486mat), dims = 2, pca = TRUE, check_duplicates = FALSE)
tsne_coords_GSE113486 <- as.data.frame(tsne_result_GSE113486$Y)
tsne_coords_GSE113486$Category <- umap_coords_GSE113486$Category

ggplot(tsne_coords_GSE113486, aes(x = V1, y = V2, color = Category)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Plot GSE113486",
       x = "t-SNE 1",
       y = "t-SNE 2") +
  scale_color_manual(values = c("cancer" = "red", "non_cancer" = "blue")) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12, face = "bold"), # Increase legend text size
    legend.title = element_text(size = 14, face = "bold") # Increase legend title size
  )

ggsave("Figures/t-sne_GSE113486.png")


GSE113486_expset <-ExpressionSet(GSE113486mat, phenoData = phenoData_GSE113486)
results_GSE113486 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                             Data = GSE113486_expset,
                                             tolerate_missed_genes = TRUE,
                                             weighted_votes = TRUE,
                                             verbose = TRUE)
knitr::kable(head(results_GSE113486))

confusion_GSE113486 <- caret::confusionMatrix(data = factor(results_GSE113486$max_score, 
                                                            levels = unique(train_object$data$Labels)),
                                              reference = factor(pData(GSE113486_expset)[,"class"], 
                                                                 levels = unique(train_object$data$Labels)),
                                              mode="everything")
print(confusion_GSE113486)


plot_confusion_matrix(confusion_GSE113486, "Confusion Matrix for GSE113486")
expr_subset_GSE113486 <- GSE113486mat[rownames(GSE11386_framed) %in% miRNAs_to_plot, ]

expr_long_GSE113486 <- expr_subset_GSE113486 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)

metaplot_GSE113486 <- as.data.frame(meta_GSE113486_h)
metaplot_GSE113486$Sample <- rownames(metaplot_GSE113486)
rownames(metaplot_GSE113486) <- NULL

predictions_GSE113486 <- data.frame(
  Sample = rownames(results_GSE113486), 
  Prediction = results_GSE113486$max_score
)
head(predictions_GSE113486)


expr_long_GSE113486 <- expr_long_GSE113486 %>%
  left_join(metaplot_GSE113486, by = "Sample") %>%
  left_join(predictions_GSE113486, by = "Sample") %>%
  rename(actual_label = class)

metaplot_GSE113486 <- metaplot_GSE113486 %>%
  rename(actual_label = class)

# Merge the actual labels with predictions
predictions_GSE113486 <- predictions_GSE113486 %>%
  left_join(metaplot_GSE113486, by = "Sample")
prediction_outcomes_GSE113486 <- predictions_GSE113486 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))
head(prediction_outcomes_GSE113486)
table(prediction_outcomes_GSE113486$Correct)


# GSE113740_Test ----------------------------------------------------------

# Fetch the GEO data
GSE113740 <- getGEO("GSE113740", GSEMatrix = TRUE)

metadata_GSE113740 = pData(phenoData(GSE113740[[1]]))
metadata_GSE113740.subset = dplyr::select(metadata_GSE113740,c(1,48))
colnames(metadata_GSE113740.subset) <- c("Title", "Description")
frequency_table_GSE113740 <- table(metadata_GSE113740.subset$Title)

# organize
metadata.subset_cancer = metadata_GSE113740.subset %>%
  dplyr::rename(state = "Description") %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))

str(GSE113740)
GSE113740_data <- GSE113740[[1]]
meta_GSE113740 <- pData(GSE113740_data)
annot_GSE113740 <- fData(GSE113740_data)

exp_GSE113740 <- exprs(GSE113740_data)
GSE113740_frame <- as.data.frame(exp_GSE113740)
rownames(GSE113740_frame) <- annot_GSE113740$miRNA_ID_LIST
GSE113740_frame <- na.omit(GSE113740_frame)

cleaned_GSE113740 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE113740_frame))
cleaned_GSE113740 <- gsub("-", "_", cleaned_GSE113740) 
cleaned_GSE113740 <- gsub(",", ".", cleaned_GSE113740) 

rownames(GSE113740_frame) <- cleaned_GSE113740
GSE113740mat <- as.matrix(GSE113740_frame)

meta_GSE113740_h <- read.csv("test/GSE113740_metadatah.csv", stringsAsFactors = TRUE, row.names = 1 )
phenoData_GSE113740 <- new("AnnotatedDataFrame",data=meta_GSE113740_h)


umap_result_GSE113740 <- umap(t(GSE113740_frame))# Plot UMAP with specific colors for cancer and non-cancer
umap_coords_GSE113740 <- as.data.frame(umap_result_GSE113740$layout)
colnames(umap_coords_GSE113740) <- c("UMAP1", "UMAP2")
umap_coords_GSE113740$State <- meta_GSE113740_h$class[match(rownames(umap_coords_GSE113740), rownames(meta_GSE113740_h))]
umap_coords_GSE113740$Category <- ifelse(umap_coords_GSE113740$State == "cancer", "cancer", "non_cancer")
unique(umap_coords_GSE113740$State)


tsne_result_GSE113740 <- Rtsne(t(GSE113740_frame), dims = 2, pca = TRUE, check_duplicates = FALSE)

# Extract t-SNE coordinates
tsne_coords_GSE113740 <- as.data.frame(tsne_result_GSE113740$Y)
tsne_coords_GSE113740$Category <- umap_coords_GSE113740$Category  # Add cancer/non-cancer labels

# Plot t-SNE
ggplot(tsne_coords_GSE113740, aes(x = V1, y = V2, color = Category)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Plot GSE113740",
       x = "t-SNE 1",
       y = "t-SNE 2") +
  scale_color_manual(values = c("cancer" = "red", "non_cancer" = "blue")) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12, face = "bold"), # Increase legend text size
    legend.title = element_text(size = 14, face = "bold") # Increase legend title size
  )

ggsave("Figures/t-sne_GSE113740.png", dpi = 600, bg = "white")


GSE113740_expset <-ExpressionSet(GSE113740mat, phenoData = phenoData_GSE113740)

results_GSE113740 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                             Data = GSE113740_expset,
                                             tolerate_missed_genes = TRUE,
                                             weighted_votes = TRUE,
                                             verbose = TRUE)
knitr::kable(head(results_GSE113740))

#write.csv(results_GSE113740, "GSE113740_resultsconf.csv")

confusion_GSE113740 <- caret::confusionMatrix(data = factor(results_GSE113740$max_score, 
                                                            levels = unique(train_object$data$Labels)),
                                              reference = factor(pData(GSE113740_expset)[,"class"], 
                                                                 levels = unique(train_object$data$Labels)),
                                              mode="everything")
print(confusion_GSE113740)


plot_confusion_matrix(confusion_GSE113740, "Confusion Matrix for GSE113740")

ggsave("Figures/confusion_GSE113740.png", dpi = 600, bg = "white")

#write.csv(as.data.frame(confusion_GSE113740$table), "final4_confusion_matrix_GSE113740.csv")
#write.csv(meta_GSE113740, "meta_all_GSE113740.csv")

expr_subset_GSE113740 <- GSE113740_frame[rownames(GSE113740_frame) %in% miRNAs_to_plot, ]

library(tidyverse)

# Convert expression data to long format for ggplot
expr_long_GSE113740 <- expr_subset_GSE113740 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)


# Create a data frame from the metadata and ensure sample names are included
metaplot_GSE113740 <- as.data.frame(meta_GSE113740_h)
metaplot_GSE113740$Sample <- rownames(metaplot_GSE113740)
rownames(metaplot_GSE113740) <- NULL  # Optionally, reset the row names

# Prepare predictions from results
predictions_GSE113740 <- data.frame(
  Sample = rownames(results_GSE113740),  # Use row names as Sample names
  Prediction = results_GSE113740$max_score  # Rename this to match your earlier context
)

# View the first few rows of the predictions data frame
head(predictions_GSE113740)


expr_long_GSE113740 <- expr_long_GSE113740 %>%
  left_join(metaplot_GSE113740, by = "Sample") %>%
  left_join(predictions_GSE113740, by = "Sample") %>%
  rename(actual_label = class)  # Ensure this matches your class label


# Rename 'class' to 'actual_label' for clarity in `metaplot_GSE113740`
metaplot_GSE113740 <- metaplot_GSE113740 %>%
  rename(actual_label = class)

# Merge the actual labels with predictions
predictions_GSE113740 <- predictions_GSE113740 %>%
  left_join(metaplot_GSE113740, by = "Sample")

# Now calculate the prediction outcome
prediction_outcomes_GSE113740 <- predictions_GSE113740 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))

# View the result to confirm
head(prediction_outcomes_GSE113740)

table(prediction_outcomes_GSE113740$Correct)


# GSE122497_test ----------------------------------------------------------
GSE122497 <- getGEO("GSE122497", GSEMatrix = TRUE)
metadata_GSE122497 = pData(phenoData(GSE122497[[1]]))
metadata_GSE122497.subset = dplyr::select(metadata_GSE122497,c(1,40))
colnames(metadata_GSE122497.subset) <- c("Title", "Description")
frequency_table_GSE122497 <- table(metadata_GSE122497.subset$Description)

# organize
metadata.subset_GSE122497 = metadata_GSE122497.subset %>%
  dplyr::rename(state = "Description") %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))
str(GSE122497)

#Expression matrix and metadata
GSE122497_data <- GSE122497[[1]]
meta_GSE122497 <- pData(GSE122497_data)
annot_GSE122497 <- fData(GSE122497_data)
exp_GSE122497 <- exprs(GSE122497_data)

GSE122497_frame <- as.data.frame(exp_GSE122497)
rownames(GSE122497_frame) <- annot_GSE122497$miRNA_ID_LIST

#For better representation
cleaned_GSE122497 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE122497_frame))
cleaned_GSE122497 <- gsub("-", "_", cleaned_GSE122497) 
cleaned_GSE122497 <- gsub(",", ".", cleaned_GSE122497) 

rownames(GSE122497_frame) <- cleaned_GSE122497
GSE122497mat <- as.matrix(GSE122497_frame)

#Make GSE122497_expression Set
meta_GSE122497_h <- read.csv("test/GSE122497_metadatah.csv", stringsAsFactors = TRUE, row.names = 1 )
phenoData_GSE122497 <- new("AnnotatedDataFrame",data=meta_GSE122497_h)
GSE122497_expset <-ExpressionSet(GSE122497mat, phenoData = phenoData_GSE122497)


# Perform t-SNE
umap_result_GSE122497 <- umap(t(GSE122497_frame))
umap_coords_GSE122497 <- as.data.frame(umap_result_GSE122497$layout)
colnames(umap_coords_GSE122497) <- c("UMAP1", "UMAP2")
umap_coords_GSE122497$State <- meta_GSE122497_h$class[match(rownames(umap_coords_GSE122497), rownames(meta_GSE122497_h))]
umap_coords_GSE122497$Category <- ifelse(umap_coords_GSE122497$State == "cancer", "cancer", "non_cancer")

tsne_result_GSE122497 <- Rtsne(t(GSE122497_frame), dims = 2, pca = TRUE, check_duplicates = FALSE)
tsne_coords_GSE122497 <- as.data.frame(tsne_result_GSE122497$Y)
tsne_coords_GSE122497$Category <- umap_coords_GSE122497$Category

# Plot t-SNE_GSE122497
ggplot(tsne_coords_GSE122497, aes(x = V1, y = V2, color = Category)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Plot GSE122497",
       x = "t-SNE 1",
       y = "t-SNE 2") +
  scale_color_manual(values = c("cancer" = "red", "non_cancer" = "blue")) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 18, face = "bold"), # Increase legend text size
    legend.title = element_text(size = 14, face = "bold") # Increase legend title size
  )

ggsave("Figures_t-sne_GSE122497.png")


results_GSE122497 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                             Data = GSE122497_expset,
                                             tolerate_missed_genes = TRUE,
                                             weighted_votes = TRUE,
                                             verbose = TRUE)
knitr::kable(head(results_GSE122497))
confusion_GSE122497 <- caret::confusionMatrix(data = factor(results_GSE122497$max_score, 
                                                            levels = unique(train_object$data$Labels)),
                                              reference = factor(pData(GSE122497_expset)[,"class"], 
                                                                 levels = unique(train_object$data$Labels)),
                                              mode="everything")
print(confusion_GSE122497)
plot_confusion_matrix(confusion_GSE122497, "Confusion Matrix for GSE122497")

expr_subset_GSE122497 <- GSE122497_frame[rownames(GSE122497_frame) %in% miRNAs_to_plot, ]
expr_long_GSE122497 <- expr_subset_GSE122497 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)

metaplot_GSE122497 <- as.data.frame(meta_GSE122497_h)
metaplot_GSE122497$Sample <- rownames(metaplot_GSE122497)
rownames(metaplot_GSE122497) <- NULL

# Prepare predictions from results
predictions_GSE122497 <- data.frame(
  Sample = rownames(results_GSE122497),  # Use row names as Sample names
  Prediction = results_GSE122497$max_score  # Rename this to match your earlier context
)
head(predictions_GSE122497)


expr_long_GSE122497 <- expr_long_GSE122497 %>%
  left_join(metaplot_GSE122497, by = "Sample") %>%
  left_join(predictions_GSE122497, by = "Sample") %>%
  rename(actual_label = class)

metaplot_GSE122497 <- metaplot_GSE122497 %>%
  rename(actual_label = class)

# Merge the actual labels with predictions
predictions_GSE122497 <- predictions_GSE122497 %>%
  left_join(metaplot_GSE122497, by = "Sample")
prediction_outcomes_GSE122497 <- predictions_GSE122497 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))
head(prediction_outcomes_GSE122497)
table(prediction_outcomes_GSE122497$Correct)


# GSE137140_Test ----------------------------------------------------------

# Fetch the GEO data
GSE137140 <- getGEO("GSE137140", GSEMatrix = TRUE)
metadata_GSE137140 = pData(phenoData(GSE137140[[1]]))
metadata_GSE137140.subset = dplyr::select(metadata_GSE137140,c(1,35))
colnames(metadata_GSE137140.subset) <- c("Title", "Description")
frequency_table_GSE137140 <- table(metadata_GSE137140.subset$Description)

# organize
metadata.subset_cancer = metadata_GSE137140.subset %>%
  dplyr::rename(state = "Description") %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))

str(GSE137140)
GSE137140_data <- GSE137140[[1]]
meta_GSE137140 <- pData(GSE137140_data)
annot_GSE137140 <- fData(GSE137140_data)
exp_GSE137140 <- exprs(GSE137140_data)
GSE137140_frame <- as.data.frame(exp_GSE137140)
rownames(GSE137140_frame) <- annot_GSE137140$miRNA_ID_LIST
GSE137140_frame <- na.omit(GSE137140_frame)

cleaned_GSE137140 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE137140_frame))
cleaned_GSE137140 <- gsub("-", "_", cleaned_GSE137140) 
cleaned_GSE137140 <- gsub(",", ".", cleaned_GSE137140) 

rownames(GSE137140_frame) <- cleaned_GSE137140
GSE137140mat <- as.matrix(GSE137140_frame)

meta_GSE137140_h <- read.csv("test/GSE137140_metadatah.csv", stringsAsFactors = TRUE, row.names = 1 )
phenoData_GSE137140 <- new("AnnotatedDataFrame",data=meta_GSE137140_h)

umap_result_GSE137140 <- umap(t(GSE137140_frame))# Plot UMAP with specific colors for cancer and non-cancer
umap_coords_GSE137140 <- as.data.frame(umap_result_GSE137140$layout)
colnames(umap_coords_GSE137140) <- c("UMAP1", "UMAP2")
umap_coords_GSE137140$State <- meta_GSE137140_h$class[match(rownames(umap_coords_GSE137140), rownames(meta_GSE137140_h))]
umap_coords_GSE137140$Category <- ifelse(umap_coords_GSE137140$State == "cancer", "cancer", "non_cancer")
unique(umap_coords_GSE137140$State)

tsne_result_GSE137140 <- Rtsne(t(GSE137140_frame), dims = 2, pca = TRUE, check_duplicates = FALSE)
tsne_coords_GSE137140 <- as.data.frame(tsne_result_GSE137140$Y)
tsne_coords_GSE137140$Category <- umap_coords_GSE137140$Category

# Plot t-SNE
ggplot(tsne_coords_GSE137140, aes(x = V1, y = V2, color = Category)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Plot GSE137140",
       x = "t-SNE 1",
       y = "t-SNE 2") +
  scale_color_manual(values = c("cancer" = "red", "non_cancer" = "blue")) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12, face = "bold"), # Increase legend text size
    legend.title = element_text(size = 14, face = "bold") # Increase legend title size
  )
ggsave("Figures/t-sne_GSE137140.png")


GSE137140_expset <-ExpressionSet(GSE137140mat, phenoData = phenoData_GSE137140)
results_GSE137140 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                             Data = GSE137140_expset,
                                             tolerate_missed_genes = TRUE,
                                             weighted_votes = TRUE,
                                             verbose = TRUE)
knitr::kable(head(results_GSE137140))
confusion_GSE137140 <- caret::confusionMatrix(data = factor(results_GSE137140$max_score, 
                                                            levels = unique(train_object$data$Labels)),
                                              reference = factor(pData(GSE137140_expset)[,"class"], 
                                                                 levels = unique(train_object$data$Labels)),
                                              mode="everything")
print(confusion_GSE137140)

plot_confusion_matrix(confusion_GSE137140, "Confusion Matrix for GSE137140")

expr_subset_GSE137140 <- GSE137140_frame[rownames(GSE137140_frame) %in% miRNAs_to_plot, ]
expr_long_GSE137140 <- expr_subset_GSE137140 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)

metaplot_GSE137140 <- as.data.frame(meta_GSE137140_h)
metaplot_GSE137140$Sample <- rownames(metaplot_GSE137140)
rownames(metaplot_GSE137140) <- NULL

# Prepare predictions from results
predictions_GSE137140 <- data.frame(
  Sample = rownames(results_GSE137140),
  Prediction = results_GSE137140$max_score
)
head(predictions_GSE137140)


expr_long_GSE137140 <- expr_long_GSE137140 %>%
  left_join(metaplot_GSE137140, by = "Sample") %>%
  left_join(predictions_GSE137140, by = "Sample") %>%
  rename(actual_label = class)  # Ensure this matches your class label
metaplot_GSE137140 <- metaplot_GSE137140 %>%
  rename(actual_label = class)

# Merge the actual labels with predictions
predictions_GSE137140 <- predictions_GSE137140 %>%
  left_join(metaplot_GSE137140, by = "Sample")
prediction_outcomes_GSE137140 <- predictions_GSE137140 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))
head(prediction_outcomes_GSE137140)


# GSE139031_Test ----------------------------------------------------------

# Fetch the GEO data
GSE139031 <- getGEO("GSE139031", GSEMatrix = TRUE)
metadata_GSE139031 = pData(phenoData(GSE139031[[1]]))
metadata_GSE139031.subset = dplyr::select(metadata_GSE139031,c(1,19))
colnames(metadata_GSE139031.subset) <- c("Title", "Description")
frequency_table_GSE139031 <- table(metadata_GSE139031.subset$Description)

# organize
metadata.subset_cancer = metadata_GSE139031.subset %>%
  dplyr::rename(state = "Description") %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))

str(GSE139031)
GSE139031_data <- GSE139031[[1]]
meta_GSE139031 <- pData(GSE139031_data)
annot_GSE139031 <- fData(GSE139031_data)

exp_GSE139031 <- exprs(GSE139031_data)
exp_139031 <- as.data.frame(exp_GSE139031)
GSE139031_frame <- as.data.frame(exp_139031)
rownames(GSE139031_frame) <- annot_GSE139031$miRNA_ID_LIST
GSE139031_frame <- na.omit(GSE139031_frame)

cleaned_GSE139031 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE139031_frame))
cleaned_GSE139031 <- gsub("-", "_", cleaned_GSE139031) 
cleaned_GSE139031 <- gsub(",", ".", cleaned_GSE139031) 
rownames(GSE139031_frame) <- cleaned_GSE139031
GSE139031mat <- as.matrix(GSE139031_frame)

meta_GSE139031_h <- read.csv("test/GSE139031_metadatah.csv", stringsAsFactors = TRUE,row.names = 1)
phenoData_GSE139031 <- new("AnnotatedDataFrame",data=meta_GSE139031_h)

umap_result_GSE139031 <- umap(t(GSE139031_frame))# Plot UMAP with specific colors for cancer and non-cancer
umap_coords_GSE139031 <- as.data.frame(umap_result_GSE139031$layout)
colnames(umap_coords_GSE139031) <- c("UMAP1", "UMAP2")

umap_coords_GSE139031$State <- meta_GSE139031_h$class[match(rownames(umap_coords_GSE139031), rownames(meta_GSE139031_h))]
umap_coords_GSE139031$Category <- ifelse(umap_coords_GSE139031$State == "cancer", "cancer", "non_cancer")
unique(umap_coords_GSE139031$State)

tsne_result_GSE139031 <- Rtsne(t(GSE139031_frame), dims = 2, pca = TRUE, check_duplicates = FALSE)
tsne_coords_GSE139031 <- as.data.frame(tsne_result_GSE139031$Y)
tsne_coords_GSE139031$Category <- umap_coords_GSE139031$Category

# Plot t-SNE
ggplot(tsne_coords_GSE139031, aes(x = V1, y = V2, color = Category)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Plot GSE139031",
       x = "t-SNE 1",
       y = "t-SNE 2") +
  scale_color_manual(values = c("cancer" = "red", "non_cancer" = "blue")) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12, face = "bold"), # Increase legend text size
    legend.title = element_text(size = 14, face = "bold") # Increase legend title size
  )

ggsave("Figures/t-sne_GSE139031.png")


GSE139031_expset <-ExpressionSet(GSE139031mat, phenoData = phenoData_GSE139031)
results_GSE139031 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                             Data = GSE139031_expset,
                                             tolerate_missed_genes = TRUE,
                                             weighted_votes = TRUE,
                                             verbose = TRUE)
knitr::kable(head(results_GSE139031))

confusion_GSE139031 <- caret::confusionMatrix(data = factor(results_GSE139031$max_score, 
                                                            levels = unique(train_object$data$Labels)),
                                              reference = factor(pData(GSE139031_expset)[,"class"], 
                                                                 levels = unique(train_object$data$Labels)),
                                              mode="everything")
print(confusion_GSE139031)
plot_confusion_matrix(confusion_GSE139031, "Confusion Matrix for GSE139031")

expr_subset_GSE139031 <- GSE139031_frame[rownames(GSE139031_frame) %in% miRNAs_to_plot, ]
expr_long_GSE139031 <- expr_subset_GSE139031 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)

metaplot_GSE139031 <- as.data.frame(meta_GSE139031_h)
metaplot_GSE139031$Sample <- rownames(metaplot_GSE139031)
rownames(metaplot_GSE139031) <- NULL  # Optionally, reset the row names

# Prepare predictions from results
predictions_GSE139031 <- data.frame(
  Sample = rownames(results_GSE139031),  # Use row names as Sample names
  Prediction = results_GSE139031$max_score  # Rename this to match your earlier context
)
head(predictions_GSE139031)


expr_long_GSE139031 <- expr_long_GSE139031 %>%
  left_join(metaplot_GSE139031, by = "Sample") %>%
  left_join(predictions_GSE139031, by = "Sample") %>%
  rename(actual_label = class)
metaplot_GSE139031 <- metaplot_GSE139031 %>%
  rename(actual_label = class)

# Merge the actual labels with predictions
predictions_GSE139031 <- predictions_GSE139031 %>%
  left_join(metaplot_GSE139031, by = "Sample")
prediction_outcomes_GSE139031 <- predictions_GSE139031 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))
head(prediction_outcomes_GSE139031)
table(prediction_outcomes_GSE139031$Correct)


# GSE164174_Test ----------------------------------------------------------

# Fetch the GEO data
GSE164174 <- getGEO("GSE164174", GSEMatrix = TRUE)
metadata_GSE164174 = pData(phenoData(GSE164174[[1]]))
metadata_GSE164174.subset = dplyr::select(metadata_GSE164174,c(1,34))
colnames(metadata_GSE164174.subset) <- c("Title", "Description")
frequency_table_GSE164174 <- table(metadata_GSE164174.subset$Description)

# organize
metadata.subset_cancer = metadata_GSE164174.subset %>%
  dplyr::rename(state = "Description") %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))

str(GSE164174)
GSE164174_data <- GSE164174[[1]]

meta_GSE164174 <- pData(GSE164174_data)
annot_GSE164174 <- fData(GSE164174_data)
exp_GSE164174 <- exprs(GSE164174_data)
frequency_table_GSE164174 <- table(metadata_GSE164174.subset$Description)

GSE164174_frame <- as.data.frame(exp_GSE164174)
rownames(GSE164174_frame) <- annot_GSE164174$miRNA_ID_LIST
GSE164174_frame <- na.omit(GSE164174_frame)

cleaned_GSE164174 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE164174_frame))
cleaned_GSE164174 <- gsub("-", "_", cleaned_GSE164174) 
cleaned_GSE164174 <- gsub(",", ".", cleaned_GSE164174) 
rownames(GSE164174_frame) <- cleaned_GSE164174
GSE164174mat <- as.matrix(GSE164174_frame)

meta_GSE164174_h <- read.csv("test/GSE164174_metadatah.csv", stringsAsFactors = TRUE, row.names = 1 )
phenoData_GSE164174 <- new("AnnotatedDataFrame",data=meta_GSE164174_h)

umap_result_GSE164174 <- umap(t(GSE164174_frame))
umap_coords_GSE164174 <- as.data.frame(umap_result_GSE164174$layout)
colnames(umap_coords_GSE164174) <- c("UMAP1", "UMAP2")
umap_coords_GSE164174$State <- meta_GSE164174_h$class[match(rownames(umap_coords_GSE164174), rownames(meta_GSE164174_h))]
umap_coords_GSE164174$Category <- ifelse(umap_coords_GSE164174$State == "cancer", "cancer", "non_cancer")
unique(umap_coords_GSE164174$State)

tsne_result_GSE164174 <- Rtsne(t(GSE164174_frame), dims = 2, pca = TRUE, check_duplicates = FALSE)
tsne_coords_GSE164174 <- as.data.frame(tsne_result_GSE164174$Y)
tsne_coords_GSE164174$Category <- umap_coords_GSE164174$Category  # Add cancer/non-cancer labels

# Plot t-SNE
ggplot(tsne_coords_GSE164174, aes(x = V1, y = V2, color = Category)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Plot GSE164174",
       x = "t-SNE 1",
       y = "t-SNE 2") +
  scale_color_manual(values = c("cancer" = "red", "non_cancer" = "blue")) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12, face = "bold"), # Increase legend text size
    legend.title = element_text(size = 14, face = "bold") # Increase legend title size
  )

ggsave("Figures/t-sne_GSE164174.png")

GSE164174_expset <-ExpressionSet(GSE164174mat, phenoData = phenoData_GSE164174)
results_GSE164174 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                             Data = GSE164174_expset,
                                             tolerate_missed_genes = TRUE,
                                             weighted_votes = TRUE,
                                             verbose = TRUE)
knitr::kable(head(results_GSE164174))

confusion_GSE164174 <- caret::confusionMatrix(data = factor(results_GSE164174$max_score, 
                                                            levels = unique(train_object$data$Labels)),
                                              reference = factor(pData(GSE164174_expset)[,"class"], 
                                                                 levels = unique(train_object$data$Labels)),
                                              mode="everything")
print(confusion_GSE164174)
plot_confusion_matrix(confusion_GSE164174, "Confusion Matrix for GSE164174")

expr_subset_GSE164174 <- GSE164174_frame[rownames(GSE164174_frame) %in% miRNAs_to_plot, ]
expr_long_GSE164174 <- expr_subset_GSE164174 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)

metaplot_GSE164174 <- as.data.frame(meta_GSE164174_h)
metaplot_GSE164174$Sample <- rownames(metaplot_GSE164174)
rownames(metaplot_GSE164174) <- NULL  # Optionally, reset the row names

# Prepare predictions from results
predictions_GSE164174 <- data.frame(
  Sample = rownames(results_GSE164174),  # Use row names as Sample names
  Prediction = results_GSE164174$max_score  # Rename this to match your earlier context
)
head(predictions_GSE164174)

expr_long_GSE164174 <- expr_long_GSE164174 %>%
  left_join(metaplot_GSE164174, by = "Sample") %>%
  left_join(predictions_GSE164174, by = "Sample") %>%
  rename(actual_label = class)  # Ensure this matches your class label

metaplot_GSE164174 <- metaplot_GSE164174 %>%
  rename(actual_label = class)

# Merge the actual labels with predictions
predictions_GSE164174 <- predictions_GSE164174 %>%
  left_join(metaplot_GSE164174, by = "Sample")
prediction_outcomes_GSE164174 <- predictions_GSE164174 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))

head(prediction_outcomes_GSE164174)
table(prediction_outcomes_GSE164174$Correct)

# GSE73002_Test -----------------------------------------------------------

# Fetch the GEO data
GSE73002 <- getGEO("GSE73002", GSEMatrix = TRUE)
metadata_GSE73002 = pData(phenoData(GSE73002[[1]]))
metadata_GSE73002.subset = dplyr::select(metadata_GSE73002,c(1,34))
colnames(metadata_GSE73002.subset) <- c("Title", "Description")
frequency_table_GSE73002 <- table(metadata_GSE73002.subset$Description)

# organize
metadata.subset_cancer = metadata_GSE73002.subset %>%
  dplyr::rename(state = "Description") %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))

str(GSE73002)
GSE73002_data <- GSE73002[[1]]

meta_GSE73002 <- pData(GSE73002_data)
annot_GSE73002 <- fData(GSE73002_data)
exp_GSE73002 <- exprs(GSE73002_data)
exp_73002 <- as.data.frame(exp_GSE73002)

rownames(exp_GSE73002) <- annot_GSE73002$miRNA_ID_LIST

#This platform has different Version of miRNA Names
miname <- annot_GSE73002$miRNA_ID_LIST

library(miRNAmeConverter)
nc = MiRNANameConverter();

# Assess the version of the miRNA names
conv <- translateMiRNAName(nc, miname, version = 21.0)

# Extract the converted miRNA names
conv_mirnas <- conv$v21.0 
mirna_rows_to_keep <- miname %in% conv_mirnas

filtered_exp_GSE73002 <- exp_GSE73002[mirna_rows_to_keep, ]

# one miRNA among 8 is missing and hence identify columns with NA values for hsa-miR-4730
na_samples_for_mirna <- is.na(filtered_exp_GSE73002["hsa-miR-4730", ])

# Get the column names (samples) with NA values
samples_with_na <- colnames(filtered_exp_GSE73002)[na_samples_for_mirna]

# Count how many NA values there are for hsa-miR-4730
na_count <- sum(na_samples_for_mirna)
print(paste("Number of NA values for hsa-miR-4730:", na_count))

# Remove only those samples (columns) where hsa-miR-4730 has NA values
filtered_exp_no_na_mirna <- filtered_exp_GSE73002[, !na_samples_for_mirna]

# Print the updated expression data for hsa-miR-4730 after removing NA samples
print(filtered_exp_no_na_mirna["hsa-miR-4730", 1:5])

filtered_exp_GSE73002 <- na.omit(filtered_exp_no_na_mirna)
metadata_filtered73002 <- metadata_GSE73002.subset %>%
  filter(!Description %in% c( "prostate disease"))
frequency_table_filtered73002 <- table(metadata_filtered73002$Description)
frequency_table_filtered73002

samples_to_keep73002 <- rownames(metadata_filtered73002)
exp_filtered_73002 <- filtered_exp_GSE73002[, colnames(filtered_exp_GSE73002) %in% samples_to_keep73002]
GSE73002_frame <- as.data.frame(exp_filtered_73002)


cleaned_GSE73002 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE73002_frame))
cleaned_GSE73002 <- gsub("-", "_", cleaned_GSE73002) 
cleaned_GSE73002 <- gsub(",", ".", cleaned_GSE73002) 
rownames(GSE73002_frame) <- cleaned_GSE73002
GSE73002mat <- as.matrix(GSE73002_frame)
GSE73002_samples <- data.frame(colnames(GSE73002_frame))

meta_GSE73002_h <- read.csv("test/GSE73002_metadatafull.csv", stringsAsFactors = TRUE,row.names = 1)
phenoData_GSE73002 <- new("AnnotatedDataFrame",data=meta_GSE73002_h)

# Extract the row names (sample IDs) from meta_GSE73002_h
sample_ids_73002 <- rownames(meta_GSE73002_h)

# Filter metadata_filtered73002 to include only rows with these sample IDs
filtered_metadata_73002 <- metadata_filtered73002[sample_ids_73002, ]

GSE73002_expset <-ExpressionSet(GSE73002mat, phenoData = phenoData_GSE73002)

umap_result_GSE73002 <- umap(t(GSE73002_frame))
umap_coords_GSE73002 <- as.data.frame(umap_result_GSE73002$layout)
colnames(umap_coords_GSE73002) <- c("UMAP1", "UMAP2")
umap_coords_GSE73002$State <- meta_GSE73002_h$class[match(rownames(umap_coords_GSE73002), rownames(meta_GSE73002_h))]
umap_coords_GSE73002$Category <- ifelse(umap_coords_GSE73002$State == "cancer", "cancer", "non_cancer")
unique(umap_coords_GSE73002$State)

tsne_result_GSE73002 <- Rtsne(t(GSE73002_frame), dims = 2, pca = TRUE, check_duplicates = FALSE)
tsne_coords_GSE73002 <- as.data.frame(tsne_result_GSE73002$Y)
tsne_coords_GSE73002$Category <- umap_coords_GSE73002$Category

# Plot t-SNE
ggplot(tsne_coords_GSE73002, aes(x = V1, y = V2, color = Category)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Plot GSE73002",
       x = "t-SNE 1",
       y = "t-SNE 2") +
  scale_color_manual(values = c("cancer" = "red", "non_cancer" = "blue")) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12, face = "bold"), # Increase legend text size
    legend.title = element_text(size = 14, face = "bold") # Increase legend title size
  )
ggsave("Figures/tsne_GSE73002_nona.png")

results_GSE73002 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                            Data = GSE73002_expset,
                                            tolerate_missed_genes = TRUE,
                                            weighted_votes = TRUE,
                                            verbose = TRUE)
knitr::kable(head(results_GSE73002))

confusion_GSE73002 <- caret::confusionMatrix(data = factor(results_GSE73002$max_score, 
                                                           levels = unique(train_object$data$Labels)),
                                             reference = factor(pData(GSE73002_expset)[,"class"], 
                                                                levels = unique(train_object$data$Labels)),
                                             mode="everything")
print(confusion_GSE73002)


plot_confusion_matrix(confusion_GSE73002, "Confusion Matrix for GSE73002")

expr_subset_GSE73002 <- GSE73002_frame[rownames(GSE73002_frame) %in% miRNAs_to_plot, ]
expr_long_GSE73002 <- expr_subset_GSE73002 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)

metaplot_GSE73002 <- as.data.frame(meta_GSE73002_h)
metaplot_GSE73002$Sample <- rownames(metaplot_GSE73002)
rownames(metaplot_GSE73002) <- NULL  # Optionally, reset the row names

predictions_GSE73002 <- data.frame(
  Sample = rownames(results_GSE73002),  # Use row names as Sample names
  Prediction = results_GSE73002$max_score  # Rename this to match your earlier context
)
head(predictions_GSE73002)


expr_long_GSE73002 <- expr_long_GSE73002 %>%
  left_join(metaplot_GSE73002, by = "Sample") %>%
  left_join(predictions_GSE73002, by = "Sample") %>%
  rename(actual_label = class)
metaplot_GSE73002 <- metaplot_GSE73002 %>%
  rename(actual_label = class)

# Merge the actual labels with predictions
predictions_GSE73002 <- predictions_GSE73002 %>%
  left_join(metaplot_GSE73002, by = "Sample")
prediction_outcomes_GSE73002 <- predictions_GSE73002 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))

head(prediction_outcomes_GSE73002)
table(prediction_outcomes_GSE73002$Correct)

# 4A_Test_Datasets_Plot_Distribution -----------------------------------------

metafull_GSE106817 <- read.csv("test/GSE106817_metadatah.csv", row.names = 1, stringsAsFactors = TRUE)
metafull_GSE112264 <- read.csv("test/GSE112264_metadatah.csv", row.names = 1, stringsAsFactors = TRUE)
metaf_GSE113486 <- read.csv("test/GSE113486_metadatah.csv", row.names = 1, stringsAsFactors = TRUE)
metafull_GSE113740 <- read.csv("test/GSE113740_metadatah.csv", row.names = 1, stringsAsFactors = TRUE)
metafull_GSE122497 <- read.csv("test/GSE122497_metadatah.csv", row.names = 1, stringsAsFactors = TRUE)
metafull_GSE137140 <- read.csv("test/GSE137140_metadatah.csv", row.names = 1, stringsAsFactors = TRUE)
metafull_GSE139031 <- read.csv("test/GSE139031_metadatah.csv", row.names = 1, stringsAsFactors = TRUE)
metafull_GSE164174 <- read.csv("test/GSE164174_metadatah.csv", row.names = 1, stringsAsFactors = TRUE)
metafull_GSE73002 <- read.csv("test/GSE73002_metadatafull.csv", row.names = 1, stringsAsFactors = TRUE)

#Duplicates_removal
metafull_GSE113486 <- metaf_GSE113486[rownames(metaf_GSE113486) %in% colnames(GSE113486mat), , drop = FALSE]

# Combine all datasets into one data frame
datasets_dist <- list(
  GSE106817 = metafull_GSE106817,
  GSE112264 = metafull_GSE112264,
  GSE113486 = metafull_GSE113486,
  GSE113740 = metafull_GSE113740,
  GSE122497 = metafull_GSE122497,
  GSE137140 = metafull_GSE137140,
  GSE139031 = metafull_GSE139031,
  GSE164174 = metafull_GSE164174,
  GSE73002 = metafull_GSE73002
)

# Add a 'dataset' column to identify each dataset
names(datasets_dist) <- c("GSE106817", "GSE112264", "GSE113486", 
                          "GSE113740", "GSE122497", "GSE137140", 
                          "GSE139031", "GSE164174", "GSE73002")

# Combine datasets into a single data frame
combined_data_dist <- do.call(rbind, lapply(names(datasets_dist), function(name) {
  data_dist <- datasets_dist[[name]]
  data_dist$dataset <- name
  return(data_dist)
}))

# Summarize the counts for each state by class and dataset
summary_data_dist <- combined_data_dist %>%
  group_by(dataset, description, class) %>%
  summarise(count = n(), .groups = "drop")

# Create the geom_tile plot
ggplot(summary_data_dist, aes(x = dataset, y = description, fill = class)) +
  geom_tile(color = "white") +
  geom_text(aes(label = count), color = "white", size = 5, fontface = "bold") +  # Bold text for counts
  scale_fill_manual(values = c("cancer" = "red", "non_cancer" = "blue"), name = "Class") +
  labs(title = "Distribution of Samples in Test Datasets",
       x = "Dataset", 
       y = "State") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black", size = 12),  # Bold x-axis labels
    axis.text.y = element_text(face = "bold", color = "black", size = 14),                        # Bold y-axis labels
    axis.title = element_text(face = "bold", color = "black", size = 16),                         # Bold axis titles
    plot.title = element_text(face = "bold", hjust = 0.5),            # Bold and centered title
    legend.title = element_text(face = "bold"),                       # Bold legend title
    legend.text = element_text(face = "bold", size = 12)                         # Bold legend text
  )

ggsave("Figures/4A_distribution_test.png", dpi = 600, bg = "white", width = 10, height = 12)


# beeswarm_test -----------------------------------------------------------

library(ggbeeswarm)
library(gridExtra)

# Define the function to create and save beeswarm plots with specific colors for labels
create_beeswarm_plots_test <- function(dataset, dataset_name) {
  # Initialize list to store beeswarm plots
  beeswarm_plots_test <- list()
  
  # Loop through each of the top 4 rules to create beeswarm plots
  for (i in seq_len(min(4, length(top_rules_train$Rule)))) {
    rule <- top_rules_train$Rule[i]
    gene1 <- top_rules_train$Gene1[i]
    gene2 <- top_rules_train$Gene2[i]
    
    # Filter data for the current rule's gene pair
    gene_data_test <- dataset %>%
      filter(Gene %in% c(gene1, gene2)) %>%
      mutate(Rule = rule)
    
    # Ensure gene_data has actual_label and Expression before plotting
    if (nrow(gene_data_test) > 0) {
      # Create the beeswarm plot with specified color for "cancer" and "no" labels
      beeswarm_plot_test <- ggplot(gene_data_test, aes(x = actual_label, y = Expression, color = actual_label)) +
        geom_beeswarm(cex = 0.5, size = 2, alpha = 0.6) +
        facet_wrap(~ Gene, scales = "free_y") +
        scale_color_manual(values = c("cancer" = "red", "non_cancer" = "blue")) +
        labs(
          title = paste("Expression Levels of miRNAs in Rule:", rule),
          x = "Actual Class",
          y = "Expression Level"
        ) +
        theme_minimal(base_size = 8) +
        theme(
          plot.title = element_text(color = "black", face = "bold", hjust = 0.5),
          axis.title.x = element_text(color = "black", face = "bold"),
          axis.title.y = element_text(color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          strip.text = element_text(color = "black", face = "bold"),
          legend.text = element_text(color = "black", face = "bold"),
          legend.title = element_text(color = "black", face = "bold"),
          legend.position = "none"
        )
      
      # Add the plot to the list
      beeswarm_plots_test[[i]] <- beeswarm_plot_test
    } else {
      message(paste("No data available for rule:", rule, "in dataset:", dataset_name))
    }
  }
  
  # Arrange the first four beeswarm plots in a 2x2 grid and save the plot
  if (length(beeswarm_plots_test) > 0) {
    grid_plot <- grid.arrange(grobs = beeswarm_plots_test, ncol = 2, nrow = 2, 
                              top = paste("Beeswarm Plots of Expression Levels for Top 4 Rules in", dataset_name))
    
    # Save the 2x2 grid of beeswarm plots as a single image in the "Figures/" directory
    ggsave(paste0("Figures/Beeswarm_Plots_", dataset_name, ".png"), grid_plot,
           dpi = 600, width = 24, height = 12, units = "cm")
  }
}


# Now apply the function to each dataset
create_beeswarm_plots_test(expr_long_GSE106817, "GSE106817")
create_beeswarm_plots_test(expr_long_GSE112264, "GSE112264")
create_beeswarm_plots_test(expr_long_GSE113486, "GSE113486")
create_beeswarm_plots_test(expr_long_GSE113740, "GSE113740")
create_beeswarm_plots_test(expr_long_GSE122497, "GSE122497")
create_beeswarm_plots_test(expr_long_GSE137140, "GSE137140")
create_beeswarm_plots_test(expr_long_GSE139031, "GSE139031")
create_beeswarm_plots_test(expr_long_GSE164174, "GSE164174")
create_beeswarm_plots_test(expr_long_GSE73002, "GSE73002")

# 4B_Confusion_Matrices_Test -------------------------------------------------

# Sample confusion matrices in list form for demonstration
confusion_matrices <- list(
  GSE122497 = confusion_GSE122497,
  GSE106817 = confusion_GSE106817,
  GSE137140 = confusion_GSE137140,
  GSE164174 = confusion_GSE164174,
  GSE113740 = confusion_GSE113740,
  GSE112264 = confusion_GSE112264,
  GSE113486 = confusion_GSE113486,
  GSE139031 = confusion_GSE139031,
  #Test_Total = confusion_test_merged)
  GSE73002 = confusion_GSE73002)

# Convert each matrix to a data frame and add a `Dataset` column
all_cm <- do.call(rbind, lapply(names(confusion_matrices), function(name) {
  cm <- as.data.frame(confusion_matrices[[name]]$table)  # assumes confusion_matrix$table format
  cm$Dataset <- name  # add a column to indicate the dataset
  cm  # return modified data frame
}))

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


ggsave("Figures/4B_Confusion_Matrix_test.png", dpi = 600, bg = "white", height = 12, width = 16)

# 4C_Test_merged -------------------------------------------------------------

# List of your datasets_m
datasets_m <- list(GSE106817 = GSE106817_expset,
                   GSE112264 = GSE112264_expset,
                   GSE113486 = GSE113486_expset,
                   GSE113740 = GSE113740_expset,
                   GSE122497 = GSE122497_expset,
                   GSE137140 = GSE137140_expset,
                   GSE139031 = GSE139031_expset,
                   GSE164174 = GSE164174_expset,
                   GSE73002 = GSE73002_expset)

# Step 1: Find common features (genes) across all datasets_m_m
common_features <- Reduce(intersect, lapply(datasets_m, function(x) featureNames(x)))

# Step 2: Subset datasets_m_m to keep only common features
datasets_m_subset <- lapply(datasets_m, function(x) x[common_features, ])

# Step 3: Merge expression matrices
expression_matrices <- lapply(datasets_m_subset, exprs)
merged_exprs <- do.call(cbind, expression_matrices)

# Step 4: Merge sample metadata and include dataset name
metadata_list <- lapply(names(datasets_m_subset), function(dataset_name) {
  meta <- pData(datasets_m_subset[[dataset_name]])
  meta$Dataset <- dataset_name  # Add dataset name as a new column
  return(meta)
})

merged_metadata <- do.call(rbind, metadata_list)

merged_meta <- merged_metadata
merged_meta$Dataset <- NULL
merged_test_mat <- as.matrix(merged_exprs)

phenoData_mergedtest <- new("AnnotatedDataFrame",data=merged_meta)
merged_expressionset <- ExpressionSet(merged_test_mat, phenoData = phenoData_mergedtest)

results_test_merged <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                               Data = merged_expressionset,
                                               tolerate_missed_genes = TRUE,
                                               weighted_votes = TRUE,
                                               verbose = TRUE)
knitr::kable(head(results_test_merged))

confusion_test_merged <- caret::confusionMatrix(data = factor(results_test_merged$max_score, 
                                                              levels = unique(train_object$data$Labels)),
                                                reference = factor(pData(merged_expressionset)[,"class"], 
                                                                   levels = unique(train_object$data$Labels)),
                                                mode="everything")
print(confusion_test_merged)

plot_confusion_matrix(confusion_test_merged, "Confusion Matrix for Test_Datasets")

ggsave("Figures/confusion_merged_Test.png")

descrip_test <- combined_data_dist %>%
  dplyr::select(class, description, dataset)

head(descrip_test)
unique(descrip_test$description)



# Extract unique descriptions
unique_descriptions <- unique(descrip_test$description)

# Define the keyword for "Benign"
benign_keyword <- "Benign"

# Initialize a color vector
colors <- character(length(unique_descriptions))

# Assign a color for "Benign" related terms
benign_color <- "#A3A500"  # Custom color for "Benign" terms

# Generate a larger number of colors for non-benign terms
num_non_benign <- sum(!sapply(unique_descriptions, function(x) grepl(benign_keyword, x, ignore.case = TRUE)))
rest_colors <- rainbow(num_non_benign)  # Generate as many colors as non-"Benign" terms

# Loop through and assign colors
color_index <- 1
for (i in seq_along(unique_descriptions)) {
  if (grepl(benign_keyword, unique_descriptions[i], ignore.case = TRUE)) {
    colors[i] <- benign_color  # Assign the same color for "Benign" terms
  } else {
    colors[i] <- rest_colors[color_index]  # Assign a unique color to non-benign terms
    color_index <- color_index + 1
  }
}

# Create a named vector for clarity
names(colors) <- unique_descriptions

# View the color assignments
colors

ref_colors <- c("cancer" = "red", "non_cancer" = "blue")

png('Figures/test_merged.png', res = 600, units = "in", width = 22, height = 10, bg = "white")

# Plot the heatmap with custom reference label colors
plot_binary_TSP(Data = merged_expressionset,               # Your data object
                classifier = classifier_train,     # Your classifier
                ref = "class",
                binary_col = c("salmon", "lightgreen", "gray"), # Binary heatmap colors
                prediction = results_test_merged,        # Your prediction data
                platform = descrip_test$description,  # Platform/study info
                platform_col = colors,     # Platform colors
                show_rule_name = TRUE,        # Show rule names
                legend = FALSE,               # Hide legend
                anno_height = 0.04,          # Annotation height
                score_height = 0.075,         # Score height
                title = "cancer_classifier_test",  # Title of the plot
                ref_col = ref_colors,
                pred_col = ref_colors,
                margin = c(0,6,0,6))
dev.off()

library(ComplexHeatmap)

platform_legend_test <- Legend(
  title = "Platform/Study",
  at = names(colors),
  legend_gp = gpar(fill = colors),   # Colors for the legend
  ncol = 2,                                  # Number of rows in the legend
  title_gp = gpar(fontsize = 14, fontface = "bold", col = "black"),  # Bold, black title
  labels_gp = gpar(fontsize = 12, fontface = "bold", col = "black")  # Bold, black labels
)

# Save the platform legend as a separate plot
png('Figures/test_merged_platform_legend.png', res = 600, units = "in", width = 14, height = 6, bg = "white")
draw(platform_legend_test)
dev.off()

# 4D_Performance_metrices_heatmap_Test ---------------------------------------

cm_metrics <- data.frame(
  Dataset = names(confusion_matrices),
  Accuracy = sapply(confusion_matrices, function(cm) {
    accuracy <- sum(diag(cm$table)) / sum(cm$table)
    return(accuracy)
  }),
  Precision = sapply(confusion_matrices, function(cm) {
    precision <- cm$table[2, 2] / sum(cm$table[2, ])
    return(precision)
  }),
  Recall = sapply(confusion_matrices, function(cm) {
    recall <- cm$table[2, 2] / sum(cm$table[, 2])
    return(recall)
  }),
  F1 = sapply(confusion_matrices, function(cm) {
    precision <- cm$table[2, 2] / sum(cm$table[2, ])
    recall <- cm$table[2, 2] / sum(cm$table[, 2])
    f1 <- 2 * (precision * recall) / (precision + recall)
    return(f1)
  }),
  Sensitivity = sapply(confusion_matrices, function(cm) {
    sensitivity <- cm$table[2, 2] / sum(cm$table[, 2])  # True Positive Rate
    return(sensitivity)
  }),
  Specificity = sapply(confusion_matrices, function(cm) {
    specificity <- cm$table[1, 1] / sum(cm$table[1, ])  # True Negative Rate
    return(specificity)
  })
)


# Load required libraries
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(scales)

# Reshape the data to long format (if not already done)
cm_metrics_long <- pivot_longer(cm_metrics, cols = -Dataset, names_to = "Metric", values_to = "Value")

# Bar plot for metrics with a custom light blue-green palette
ggplot(cm_metrics_long, aes(x = Dataset, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Performance Metrics for Multiple Datasets",
       x = "Dataset", 
       y = "Metric Value") +
  scale_fill_manual(values = c("#AEDFF7", "#7FB3D5", "#2E86C1", "#85C1E9", "#73C6B6", "#1ABC9C")) +  # Custom light color palette
  theme_minimal() +
  theme(
    legend.position = "top",  # Move legend to the top
    legend.title = element_blank(),  # Remove the legend title
    axis.text = element_text(size = 10),  # Improve axis text readability
    strip.text = element_text(size = 12, face = "bold"),  # Improve facet text size
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Title size and alignment
    axis.title = element_text(size = 12),  # Axis title size
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text for better readability
  ) +
  scale_y_continuous(labels = percent)  # Convert y-axis to percentage representation

# Heatmap plot with a lighter blue gradient for readability
ggplot(cm_metrics_long, aes(x = Metric, y = Dataset, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = percent(Value, accuracy = 0.1)), color = "black", size = 3) +  # Add percentage labels
  scale_fill_gradient(low = "#E0F7FA", high = "#0288D1") +  # Custom gradient from light blue to darker blue
  labs(title = "Rules-based classifier Metrics - Test Datasets", 
       x = "Metric", 
       y = "Dataset", 
       fill = "Metric Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Center the title
        axis.title = element_text(size = 12),  # Improve axis title size
        strip.text = element_text(size = 12, face = "bold")) +  # Bold and larger facet labels
  scale_x_discrete(expand = c(0.1, 0.1))  # Add some spacing for x-axis labels


library(scales)  # For percent() function

ggplot(cm_metrics_long, aes(x = Metric, y = Dataset, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = percent(Value, accuracy = 0.1)), 
            color = "black", size = 5, fontface = "bold") +  # Bold percentage labels
  scale_fill_gradient(low = "#E0F7FA", high = "#0288D1") +  # Custom gradient from light blue to darker blue
  labs(title = "Rules-based Classifier Metrics - Test Datasets", 
       x = "Metric", 
       y = "Dataset", 
       fill = "Metric Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10,color = "black"),  # Bold x-axis labels
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),                         # Bold y-axis labels
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),             # Bold and centered title
    axis.title.x = element_text(size = 12, face = "bold", color = "black"),                        # Bold x-axis title
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),                        # Bold y-axis title
    strip.text = element_text(size = 12, face = "bold"),                          # Bold facet labels
    legend.text = element_text(size = 10, face = "bold"),                         # Bold legend text
    legend.title = element_text(size = 12, face = "bold")                         # Bold legend title
  ) +
  scale_x_discrete(expand = c(0.1, 0.1))  # Add some spacing for x-axis labels

ggsave("heatmap_metrices_bold_ordered.png", dpi = 600, bg = "white", height = 6, width = 8)

# Order of the levels for the Dataset column
cm_metrics_long$Dataset <- factor(cm_metrics_long$Dataset, 
                                  levels = rev(c("GSE106817", "GSE112264", "GSE113486", 
                                                 "GSE113740", "GSE122497", "GSE137140", 
                                                 "GSE139031", "GSE164174", "GSE73002")))

# Verify levels (optional step to debug)
print(levels(cm_metrics_long$Dataset))

# Create the plot
ggplot(cm_metrics_long, aes(x = Metric, y = Dataset, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = percent(Value, accuracy = 0.1)), 
            color = "black", size = 3, fontface = "bold") +  # Bold percentage labels
  scale_fill_gradient(low = "#E0F7FA", high = "#0288D1") +  # Custom gradient from light blue to darker blue
  labs(title = "Rules-based Classifier Metrics - Test Datasets", 
       x = "Metric", 
       y = "Dataset", 
       fill = "Metric Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),  # Bold x-axis labels
    axis.text.y = element_text(face = "bold", size = 10),                         # Bold y-axis labels
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),             # Bold and centered title
    axis.title.x = element_text(size = 12, face = "bold"),                        # Bold x-axis title
    axis.title.y = element_text(size = 12, face = "bold"),                        # Bold y-axis title
    strip.text = element_text(size = 12, face = "bold"),                          # Bold facet labels
    legend.text = element_text(size = 10, face = "bold"),                         # Bold legend text
    legend.title = element_text(size = 12, face = "bold")                         # Bold legend title
  ) +
  scale_x_discrete(expand = c(0.1, 0.1))  # Add some spacing for x-axis labels

library(scales)  # For percent() function

ggplot(cm_metrics_long, aes(x = Metric, y = Dataset, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = percent(Value, accuracy = 0.1)), 
            color = "black", size = 5, fontface = "bold") +  # Bold percentage labels
  scale_fill_gradient(low = "lavender", high = "lightblue") +  # Custom gradient from light blue to darker blue
  labs(title = "Rules-based Classifier Metrics - Test Datasets", 
       x = "Metric", 
       y = "Dataset", 
       fill = "Metric") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 16, color = "black"),  # Bold x-axis labels
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


ggsave("Figures/4C_wb_heatmap_metrices_bold_short.png")


# 4E_ROC_Test ----------------------------------------------------------------

true_labels_122497 <- as.numeric(factor(meta_GSE122497_h$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_122497 <- results_GSE122497$cancer
roc_GSE122497 <- roc(true_labels_122497, predicted_probs_122497)
plot(roc_GSE122497, col = "blue", main = "ROC Curve for GSE122497", lwd = 2)
auc_value_122497 <- auc(roc_GSE122497)
text(0.6, 0.2, paste("AUC =", round(auc_value_122497, 4)), col = "blue", cex = 1.2)

true_labels_106817 <- as.numeric(factor(meta_GSE106817_h$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_106817 <- results_GSE106817$cancer
roc_GSE106817 <- roc(true_labels_106817, predicted_probs_106817)
plot(roc_GSE106817, col = "blue", main = "ROC Curve for GSE106817", lwd = 2)
auc_value_106817 <- auc(roc_GSE106817)
text(0.6, 0.2, paste("AUC =", round(auc_value_106817, 4)), col = "blue", cex = 1.2)


true_labels_137140 <- as.numeric(factor(meta_GSE137140_h$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_137140 <- results_GSE137140$cancer
roc_GSE137140 <- roc(true_labels_137140, predicted_probs_137140)
plot(roc_GSE137140, col = "blue", main = "ROC Curve for GSE137140", lwd = 2)
auc_value_137140 <- auc(roc_GSE137140)
text(0.6, 0.2, paste("AUC =", round(auc_value_137140, 4)), col = "blue", cex = 1.2)

true_labels_164174 <- as.numeric(factor(meta_GSE164174_h$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_164174 <- results_GSE164174$cancer
roc_GSE164174 <- roc(true_labels_164174, predicted_probs_164174)
plot(roc_GSE164174, col = "blue", main = "ROC Curve for GSE164174", lwd = 2)
auc_value_164174 <- auc(roc_GSE164174)
text(0.6, 0.2, paste("AUC =", round(auc_value_164174, 4)), col = "blue", cex = 1.2)

true_labels_113740 <- as.numeric(factor(meta_GSE113740_h$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_113740 <- results_GSE113740$cancer
roc_GSE113740 <- roc(true_labels_113740, predicted_probs_113740)
plot(roc_GSE113740, col = "blue", main = "ROC Curve for GSE113740", lwd = 2)
auc_value_113740 <- auc(roc_GSE113740)
text(0.6, 0.2, paste("AUC =", round(auc_value_113740, 4)), col = "blue", cex = 1.2)


true_labels_112264 <- as.numeric(factor(meta_GSE112264_h$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_112264 <- results_GSE112264$cancer
roc_GSE112264 <- roc(true_labels_112264, predicted_probs_112264)
plot(roc_GSE112264, col = "blue", main = "ROC Curve for GSE112264", lwd = 2)
auc_value_112264 <- auc(roc_GSE112264)
text(0.6, 0.2, paste("AUC =", round(auc_value_112264, 4)), col = "blue", cex = 1.2)


true_labels_113486 <- as.numeric(factor(meta_GSE113486_h$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_113486 <- results_GSE113486$cancer
roc_GSE113486 <- roc(true_labels_113486, predicted_probs_113486)
plot(roc_GSE113486, col = "blue", main = "ROC Curve for GSE113486", lwd = 2)
auc_value_113486 <- auc(roc_GSE113486)
text(0.6, 0.2, paste("AUC =", round(auc_value_113486, 4)), col = "blue", cex = 1.2)

true_labels_139031 <- as.numeric(factor(meta_GSE139031_h$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_139031 <- results_GSE139031$cancer
roc_GSE139031 <- roc(true_labels_139031, predicted_probs_139031)
plot(roc_GSE139031, col = "blue", main = "ROC Curve for GSE139031", lwd = 2)
auc_value_139031 <- auc(roc_GSE139031)
text(0.6, 0.2, paste("AUC =", round(auc_value_139031, 4)), col = "blue", cex = 1.2)


#true_labels_Test_Total <- as.numeric(factor(merged_metadata$class, levels = c("non_cancer", "cancer"))) - 1
#predicted_probs_Test_Total <- results_test_merged$cancer
#roc_Test_Total <- roc(true_labels_Test_Total, predicted_probs_Test_Total)
#plot(roc_Test_Total, col = "blue", main = "ROC Curve for Test_Total", lwd = 2)
#auc_value_Test_Total <- auc(roc_Test_Total)
#text(0.6, 0.2, paste("AUC =", round(auc_value_Test_Total, 4)), col = "blue", cex = 1.2)


true_labels_73002 <- as.numeric(factor(meta_GSE73002_h$class, levels = c("non_cancer", "cancer"))) - 1
predicted_probs_73002 <- results_GSE73002$cancer
roc_GSE73002 <- roc(true_labels_73002, predicted_probs_73002)
plot(roc_GSE73002, col = "blue", main = "ROC Curve for GSE73002", lwd = 2)
auc_value_73002 <- auc(roc_GSE73002)
text(0.6, 0.2, paste("AUC =", round(auc_value_73002, 4)), col = "blue", cex = 1.2)


# Step 1: Prepare the ROC data as you did before
roc_data <- list(
  GSE122497 = roc_GSE122497,
  GSE106817 = roc_GSE106817,
  GSE137140 = roc_GSE137140,
  GSE164174 = roc_GSE164174,
  GSE113740 = roc_GSE113740,
  GSE112264 = roc_GSE112264,
  GSE113486 = roc_GSE113486,
  GSE139031 = roc_GSE139031,
  #Test_Total = roc_Test_Total
  GSE73002 = roc_GSE73002
)

# Convert each ROC object to a data frame and add a `Dataset` column
all_roc <- do.call(rbind, lapply(names(roc_data), function(name) {
  roc <- roc_data[[name]]
  df <- data.frame(
    FPR = 1 - roc$specificities,  # False Positive Rate
    TPR = roc$sensitivities,      # True Positive Rate
    Threshold = roc$thresholds    # Thresholds for ROC
  )
  df$Dataset <- name  # add dataset name as a column
  
  # Calculate AUC for each ROC object
  auc_value <- auc(roc)  # Compute the AUC for the current ROC curve
  
  df$AUC <- auc_value  # Add the AUC value to the dataframe for later use
  return(df)
}))

# Step 2: Plot ROC Curves with ggplot
library(ggplot2)

ggplot(all_roc, aes(x = FPR, y = TPR, color = Dataset)) +
  geom_line() +
  
  # Add AUC value near the top-right corner (adjusted for visibility)
  geom_label(data = subset(all_roc, FPR == max(FPR) & TPR == max(TPR)), 
             aes(label = paste("AUC = ", format(AUC, digits = 4))),
             color = "white", size = 2.5, fontface = "bold", 
             fill = "black", label.padding = unit(0.5, "lines"),
             label.r = unit(0.15, "lines"), hjust = 1.5, vjust = 1.5) +
  
  # Add Dataset name in large font at the top-left corner
  geom_text(data = subset(all_roc, FPR == min(FPR) & TPR == max(TPR)),
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
auc_positions <- all_roc %>%
  group_by(Dataset) %>%
  summarise(
    FPR_AUC = max(FPR), TPR_AUC = max(TPR),  # Position for AUC label
    AUC_Value = unique(AUC)  # Extract the AUC value for the dataset
  ) %>%
  arrange(Dataset)  # Ensure the datasets are in order

# Ensure Dataset is a factor with the desired order
all_roc$Dataset <- factor(all_roc$Dataset, levels = auc_positions$Dataset)

# Plot
ggplot(all_roc, aes(x = FPR, y = TPR, color = Dataset)) +
  geom_line(size = 1) +
  
  # Add AUC value near the top-right corner for each dataset
  geom_label(data = auc_positions,
             aes(x = FPR_AUC, y = TPR_AUC,
                 label = paste("AUC = ", format(AUC_Value, digits = 4))),
             inherit.aes = FALSE,  # Prevent aesthetics inheritance
             color = "white", fill = "black", size = 3, fontface = "bold",
             label.padding = unit(0.5, "lines"), hjust = 1, vjust = 1) +
  
  theme_minimal() +
  labs(title = "ROC Curves for Test Datasets",
       x = "False Positive Rate", y = "True Positive Rate") +
  scale_color_viridis_d() +
  facet_wrap(~ Dataset, ncol = 3) +  # Create a 3-column grid layout
  theme(
    legend.position = "none",  # Remove legend as we're using color for dataset
    strip.text = element_text(size = 12, face = "bold"),  # Bold facet labels
    axis.title.x = element_text(size = 14, face = "bold"),  # Bold x-axis title
    axis.title.y = element_text(size = 14, face = "bold"),  # Bold y-axis title
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),  # Bold x-axis text
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),  # Bold y-axis text
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Bold and centered title
    plot.margin = margin(5, 5, 5, 5)  # Adjust plot margins
  )

ggsave("Figures/4D_roc_plots_bold_benign.png", dpi = 600, bg = "white", height = 8, width = 8)
