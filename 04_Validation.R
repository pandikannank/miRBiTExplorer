
# 04_Validation -----------------------------------------------------------


#This code will reproduce the pre-processing of Validation Data and the classifier performance
#Fig-5A,5B,5D,5F,5H can be generated.

# Validation_Datasets -----------------------------------------------------

# GSE110317_Validation ----------------------------------------------------

# Fetch the GEO data
GSE110317 <- getGEO("GSE110317", GSEMatrix = TRUE)
metadata_GSE110317 = pData(phenoData(GSE110317[[1]]))
metadata_GSE110317.subset = dplyr::select(metadata_GSE110317,c(1,39))
colnames(metadata_GSE110317.subset) <- c("Title", "Description")
frequency_table_GSE110317 <- table(metadata_GSE110317.subset$Description)

# organize
metadata.subset_cancer_GSE110317 = metadata_GSE110317.subset %>%
  dplyr::rename(state = "Description") %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))

str(GSE110317)

GSE110317_data <- GSE110317[[1]]
meta_GSE110317 <- pData(GSE110317_data)
annot_GSE110317 <- fData(GSE110317_data)
exp_GSE110317 <- exprs(GSE110317_data)
rownames(exp_GSE110317) <- annot_GSE110317$miRNA_ID_LIST

GSE110317_frame <- as.data.frame(exp_GSE110317)
GSE110317_frame <- na.omit(GSE110317_frame)

cleaned_GSE110317 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE110317_frame))
cleaned_GSE110317 <- gsub("-", "_", cleaned_GSE110317) 
cleaned_GSE110317 <- gsub(",", ".", cleaned_GSE110317) 
rownames(GSE110317_frame) <- cleaned_GSE110317
GSE110317mat <- as.matrix(GSE110317_frame)



meta_GSE110317_h <- read.csv("validation/GSE110317_metatest.csv", stringsAsFactors = TRUE,row.names = 1)
phenoData_GSE110317 <- new("AnnotatedDataFrame",data=meta_GSE110317_h)

GSE110317_expset <-ExpressionSet(GSE110317mat, phenoData = phenoData_GSE110317)


results_GSE110317 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                             Data = GSE110317_expset,
                                             tolerate_missed_genes = TRUE,
                                             weighted_votes = TRUE,
                                             verbose = TRUE)
knitr::kable(head(results_GSE110317))

confusion_GSE110317 <- caret::confusionMatrix(data = factor(results_GSE110317$max_score, 
                                                            levels = unique(train_object$data$Labels)),
                                              reference = factor(pData(GSE110317_expset)[,"class"], 
                                                                 levels = unique(train_object$data$Labels)),
                                              mode="everything")
print(confusion_GSE110317)
plot_confusion_matrix(confusion_GSE110317, "Confusion Matrix for GSE110317")


expr_subset_GSE110317 <- GSE110317_frame[rownames(GSE110317_frame) %in% miRNAs_to_plot, ]
expr_long_GSE110317 <- expr_subset_GSE110317 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)

metaplot_GSE110317 <- as.data.frame(meta_GSE110317_h)
metaplot_GSE110317$Sample <- rownames(metaplot_GSE110317)
rownames(metaplot_GSE110317) <- NULL

# Prepare predictions from results
predictions_GSE110317 <- data.frame(
  Sample = rownames(results_GSE110317), 
  Prediction = results_GSE110317$max_score
)
head(predictions_GSE110317)


expr_long_GSE110317 <- expr_long_GSE110317 %>%
  left_join(metaplot_GSE110317, by = "Sample") %>%
  left_join(predictions_GSE110317, by = "Sample") %>%
  rename(actual_label = class)
metaplot_GSE110317 <- metaplot_GSE110317 %>%
  rename(actual_label = class)

# Merge the actual labels with predictions
predictions_GSE110317 <- predictions_GSE110317 %>%
  left_join(metaplot_GSE110317, by = "Sample")
prediction_outcomes_GSE110317 <- predictions_GSE110317 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))

head(prediction_outcomes_GSE110317)
table(prediction_outcomes_GSE110317$Correct)



# GSE110651_Validation ----------------------------------------------------

# Fetch the GEO data
GSE110651 <- getGEO("GSE110651", GSEMatrix = TRUE)
metadata_GSE110651 = pData(phenoData(GSE110651[[1]]))
metadata_GSE110651.subset = dplyr::select(metadata_GSE110651,c(2,35))
colnames(metadata_GSE110651.subset) <- c("Title", "Description")
frequency_table_GSE110651 <- table(metadata_GSE110651.subset$Description)

# organize
metadata.subset_cancer_GSE110651 = metadata_GSE110651.subset %>%
  dplyr::rename(state = "Description") %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))

str(GSE110651)
GSE110651_data <- GSE110651[[1]]
meta_GSE110651 <- pData(GSE110651_data)
annot_GSE110651 <- fData(GSE110651_data)
exp_GSE110651 <- exprs(GSE110651_data)
rownames(exp_GSE110651) <- annot_GSE110651$miRNA_ID_LIST
GSE110651_frame <- as.data.frame(exp_GSE110651)
GSE110651_frame <- na.omit(GSE110651_frame)

cleaned_GSE110651 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE110651_frame))
cleaned_GSE110651 <- gsub("-", "_", cleaned_GSE110651) 
cleaned_GSE110651 <- gsub(",", ".", cleaned_GSE110651) 
rownames(GSE110651_frame) <- cleaned_GSE110651
GSE110651mat <- as.matrix(GSE110651_frame)

meta_GSE110651_h <- read.csv("validation/GSE110651_metatest.csv", stringsAsFactors = TRUE,row.names = 1)
phenoData_GSE110651 <- new("AnnotatedDataFrame",data=meta_GSE110651_h)
GSE110651_expset <-ExpressionSet(GSE110651mat, phenoData = phenoData_GSE110651)


results_GSE110651 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                             Data = GSE110651_expset,
                                             tolerate_missed_genes = TRUE,
                                             weighted_votes = TRUE,
                                             verbose = TRUE)
knitr::kable(head(results_GSE110651))
confusion_GSE110651 <- caret::confusionMatrix(data = factor(results_GSE110651$max_score, 
                                                            levels = unique(train_object$data$Labels)),
                                              reference = factor(pData(GSE110651_expset)[,"class"], 
                                                                 levels = unique(train_object$data$Labels)),
                                              mode="everything")
print(confusion_GSE110651)
plot_confusion_matrix(confusion_GSE110651, "Confusion Matrix for GSE110651")

expr_subset_GSE110651 <- GSE110651_frame[rownames(GSE110651_frame) %in% miRNAs_to_plot, ]
expr_long_GSE110651 <- expr_subset_GSE110651 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)

metaplot_GSE110651 <- as.data.frame(meta_GSE110651_h)
metaplot_GSE110651$Sample <- rownames(metaplot_GSE110651)
rownames(metaplot_GSE110651) <- NULL

# Prepare predictions from results
predictions_GSE110651 <- data.frame(
  Sample = rownames(results_GSE110651),
  Prediction = results_GSE110651$max_score
)
head(predictions_GSE110651)


expr_long_GSE110651 <- expr_long_GSE110651 %>%
  left_join(metaplot_GSE110651, by = "Sample") %>%
  left_join(predictions_GSE110651, by = "Sample") %>%
  rename(actual_label = class)
metaplot_GSE110651 <- metaplot_GSE110651 %>%
  rename(actual_label = class)

predictions_GSE110651 <- predictions_GSE110651 %>%
  left_join(metaplot_GSE110651, by = "Sample")
prediction_outcomes_GSE110651 <- predictions_GSE110651 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))

head(prediction_outcomes_GSE110651)
table(prediction_outcomes_GSE110651$Correct)


# GSE117064_Validation ----------------------------------------------------

# Fetch the GEO data
GSE117064 <- getGEO("GSE117064", GSEMatrix = TRUE)
metadata_GSE117064 = pData(phenoData(GSE117064[[1]]))
metadata_GSE117064.subset = dplyr::select(metadata_GSE117064,c(1,8))
colnames(metadata_GSE117064.subset) <- c("Title", "Description")
frequency_table_GSE117064 <- table(metadata_GSE117064.subset$Description)

# organize
metadata.subset_cancer_GSE117064 = metadata_GSE117064.subset %>%
  dplyr::rename(state = "Description") %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))

str(GSE117064)
GSE117064_data <- GSE117064[[1]]
meta_GSE117064 <- pData(GSE117064_data)
annot_GSE117064 <- fData(GSE117064_data)
exp_GSE117064 <- exprs(GSE117064_data)
rownames(exp_GSE117064) <- annot_GSE117064$miRNA_ID_LIST
GSE117064_frame <- as.data.frame(exp_GSE117064)
GSE117064_frame <- na.omit(GSE117064_frame)

cleaned_GSE117064 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE117064_frame))
cleaned_GSE117064 <- gsub("-", "_", cleaned_GSE117064) 
cleaned_GSE117064 <- gsub(",", ".", cleaned_GSE117064) 

rownames(GSE117064_frame) <- cleaned_GSE117064
GSE117064mat <- as.matrix(GSE117064_frame)

meta_GSE117064_h <- read.csv("validation/GSE117064_metatest.csv", stringsAsFactors = TRUE, row.names = 1)
phenoData_GSE117064 <- new("AnnotatedDataFrame",data=meta_GSE117064_h)

GSE117064_expset <-ExpressionSet(GSE117064mat, phenoData = phenoData_GSE117064)

#Duplicates present in this dataset
# Determine the maximum length
max_length <- max(length(colnames(GSE117064_frame)), length(colnames(GSE113740_frame)))

# Calculate the lengths of the column names
len_GSE117064 <- length(colnames(GSE117064_frame))
len_GSE113740 <- length(colnames(GSE113740_frame))

# Determine the maximum length
max_length <- max(len_GSE117064, len_GSE113740)

# Pad the columns with NA to match the maximum length
GSE117064_padded <- c(colnames(GSE117064_frame), rep(NA, max_length - len_GSE117064))
GSE113740_padded <- c(colnames(GSE113740_frame), rep(NA, max_length - len_GSE113740))

# Create the data frame
sampleNames <- data.frame(
  GSE117064 = GSE117064_padded,
  GSE113740 = GSE113740_padded
)

# View the data frame
print(sampleNames)

common_samples <- intersect(sampleNames$GSE117064, sampleNames$GSE113740)
print(common_samples)

if (!is.data.frame(meta_GSE117064_h)) {
  meta_GSE117064_h <- as.data.frame(meta_GSE117064_h)
}

meta_GSE117064_h_filtered <- meta_GSE117064_h[!(rownames(meta_GSE117064_h) %in% common_samples), , drop = FALSE]
GSE117064mat_filtered <- GSE117064mat[, !(colnames(GSE117064mat) %in% common_samples)]

meta_GSE117064_h_filtered_df <- data.frame(
  sampleNames = rownames(meta_GSE117064_h_filtered),  # Add sample names as a column
  class = as.character(meta_GSE117064_h_filtered[, 1]) # Convert factor to character if needed
)
rownames(meta_GSE117064_h_filtered_df) <- meta_GSE117064_h_filtered_df$sampleNames  # Set sample names as rownames
meta_GSE117064_h <- meta_GSE117064_h_filtered_df 
meta_GSE117064_h$sampleNames <- NULL

phenoData_GSE117064_filtered <- new("AnnotatedDataFrame", data = meta_GSE117064_h)

GSE117064_expset_filtered <- ExpressionSet(
  assayData = GSE117064mat_filtered,
  phenoData = phenoData_GSE117064_filtered
)

metadata_GSE117064.subset_filt <- metadata_GSE117064.subset[!rownames(metadata_GSE117064.subset) %in% common_samples, ]
head(metadata_GSE117064.subset_filt)

results_GSE117064 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                             Data = GSE117064_expset_filtered,
                                             tolerate_missed_genes = TRUE,
                                             weighted_votes = TRUE,
                                             verbose = TRUE)
knitr::kable(head(results_GSE117064))

confusion_GSE117064 <- caret::confusionMatrix(data = factor(results_GSE117064$max_score, 
                                                            levels = unique(train_object$data$Labels)),
                                              reference = factor(pData(GSE117064_expset_filtered)[,"class"], 
                                                                 levels = unique(train_object$data$Labels)),
                                              mode="everything")
print(confusion_GSE117064)

plot_confusion_matrix(confusion_GSE117064, "Confusion Matrix for GSE117064")

GSE117064_frame <- data.frame(GSE117064mat_filtered)
expr_subset_GSE117064 <- GSE117064_frame[rownames(GSE117064_frame) %in% miRNAs_to_plot, ]

expr_long_GSE117064 <- expr_subset_GSE117064 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)

metaplot_GSE117064 <- as.data.frame(meta_GSE117064_h)
metaplot_GSE117064$Sample <- rownames(metaplot_GSE117064)
rownames(metaplot_GSE117064) <- NULL

# Prepare predictions from results
predictions_GSE117064 <- data.frame(
  Sample = rownames(results_GSE117064),
  Prediction = results_GSE117064$max_score
)
head(predictions_GSE117064)


expr_long_GSE117064 <- expr_long_GSE117064 %>%
  left_join(metaplot_GSE117064, by = "Sample") %>%
  left_join(predictions_GSE117064, by = "Sample") %>%
  rename(actual_label = class)

metaplot_GSE117064 <- metaplot_GSE117064 %>%
  rename(actual_label = class)

# Merge the actual labels with predictions
predictions_GSE117064 <- predictions_GSE117064 %>%
  left_join(metaplot_GSE117064, by = "Sample")

prediction_outcomes_GSE117064 <- predictions_GSE117064 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))

head(prediction_outcomes_GSE117064)
table(prediction_outcomes_GSE117064$Correct)


# GSE120584_Validation ----------------------------------------------------

# Fetch the GEO data
GSE120584 <- getGEO("GSE120584", GSEMatrix = TRUE)

metadata_GSE120584 = pData(phenoData(GSE120584[[1]]))
metadata_GSE120584.subset = dplyr::select(metadata_GSE120584,c(1,35))
colnames(metadata_GSE120584.subset) <- c("Title", "Description")
frequency_table_GSE120584 <- table(metadata_GSE120584.subset$Description)

# organize
metadata.subset_cancer_GSE120584 = metadata_GSE120584.subset %>%
  dplyr::rename(state = "Description") %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))

str(GSE120584)
GSE120584_data <- GSE120584[[1]]
meta_GSE120584 <- pData(GSE120584_data)
annot_GSE120584 <- fData(GSE120584_data)
exp_GSE120584 <- exprs(GSE120584_data)
rownames(exp_GSE120584) <- annot_GSE120584$miRNA_ID_LIST

GSE120584_frame <- as.data.frame(exp_GSE120584)
GSE120584_frame <- na.omit(GSE120584_frame)
cleaned_GSE120584 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE120584_frame))
cleaned_GSE120584 <- gsub("-", "_", cleaned_GSE120584) 
cleaned_GSE120584 <- gsub(",", ".", cleaned_GSE120584) 
rownames(GSE120584_frame) <- cleaned_GSE120584
GSE120584mat <- as.matrix(GSE120584_frame)

meta_GSE120584_h <- read.csv("validation/GSE120584_metadatatest.csv", stringsAsFactors = TRUE,row.names = 1)
phenoData_GSE120584 <- new("AnnotatedDataFrame",data=meta_GSE120584_h)

GSE120584_expset <-ExpressionSet(GSE120584mat, phenoData = phenoData_GSE120584)
results_GSE120584 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                             Data = GSE120584_expset,
                                             tolerate_missed_genes = TRUE,
                                             weighted_votes = TRUE,
                                             verbose = TRUE)
knitr::kable(head(results_GSE120584))

confusion_GSE120584 <- caret::confusionMatrix(data = factor(results_GSE120584$max_score, 
                                                            levels = unique(train_object$data$Labels)),
                                              reference = factor(pData(GSE120584_expset)[,"class"], 
                                                                 levels = unique(train_object$data$Labels)),
                                              mode="everything")
print(confusion_GSE120584)
plot_confusion_matrix(confusion_GSE120584, "Confusion Matrix for GSE120584")

expr_subset_GSE120584 <- GSE120584_frame[rownames(GSE120584_frame) %in% miRNAs_to_plot, ]
expr_long_GSE120584 <- expr_subset_GSE120584 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)

metaplot_GSE120584 <- as.data.frame(meta_GSE120584_h)
metaplot_GSE120584$Sample <- rownames(metaplot_GSE120584)
rownames(metaplot_GSE120584) <- NULL

# Prepare predictions from results
predictions_GSE120584 <- data.frame(
  Sample = rownames(results_GSE120584),
  Prediction = results_GSE120584$max_score
)
head(predictions_GSE120584)

expr_long_GSE120584 <- expr_long_GSE120584 %>%
  left_join(metaplot_GSE120584, by = "Sample") %>%
  left_join(predictions_GSE120584, by = "Sample") %>%
  rename(actual_label = class)

metaplot_GSE120584 <- metaplot_GSE120584 %>%
  rename(actual_label = class)

predictions_GSE120584 <- predictions_GSE120584 %>%
  left_join(metaplot_GSE120584, by = "Sample")

prediction_outcomes_GSE120584 <- predictions_GSE120584 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))

head(prediction_outcomes_GSE120584)
table(prediction_outcomes_GSE120584$Correct)

# GSE134108_Validation ----------------------------------------------------

# Fetch the GEO data
GSE134108 <- getGEO("GSE134108", GSEMatrix = TRUE)
metadata_GSE134108 = pData(phenoData(GSE134108[[2]]))
metadata_GSE134108.subset = dplyr::select(metadata_GSE134108,c(2,40))
colnames(metadata_GSE134108.subset) <- c("Title", "Description")
frequency_table_GSE134108 <- table(metadata_GSE134108.subset$Description)

# organize
metadata.subset_cancer_GSE134108 = metadata_GSE134108.subset %>%
  dplyr::rename(state = "Description") %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))

str(GSE134108)
GSE134108_data <- GSE134108[[2]]
meta_GSE134108 <- pData(GSE134108_data)
annot_GSE134108 <- fData(GSE134108_data)
exp_GSE134108 <- exprs(GSE134108_data)
rownames(exp_GSE134108) <- annot_GSE134108$miRNA_ID_LIST

GSE134108_frame <- as.data.frame(exp_GSE134108)
GSE134108_frame <- na.omit(GSE134108_frame)

cleaned_GSE134108 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE134108_frame))
cleaned_GSE134108 <- gsub("-", "_", cleaned_GSE134108) 
cleaned_GSE134108 <- gsub(",", ".", cleaned_GSE134108) 
rownames(GSE134108_frame) <- cleaned_GSE134108
GSE134108mat <- as.matrix(GSE134108_frame)

meta_GSE134108_h <- read.csv("validation/GSE134108_metatest.csv", stringsAsFactors = TRUE,row.names = 1)
phenoData_GSE134108 <- new("AnnotatedDataFrame",data=meta_GSE134108_h)

GSE134108_expset <-ExpressionSet(GSE134108mat, phenoData = phenoData_GSE134108)
results_GSE134108 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                             Data = GSE134108_expset,
                                             tolerate_missed_genes = TRUE,
                                             weighted_votes = TRUE,
                                             verbose = TRUE)
knitr::kable(head(results_GSE134108))

confusion_GSE134108 <- caret::confusionMatrix(data = factor(results_GSE134108$max_score, 
                                                            levels = unique(train_object$data$Labels)),
                                              reference = factor(pData(GSE134108_expset)[,"class"], 
                                                                 levels = unique(train_object$data$Labels)),
                                              mode="everything")
print(confusion_GSE134108)
plot_confusion_matrix(confusion_GSE134108, "Confusion Matrix for GSE134108")

expr_subset_GSE134108 <- GSE134108_frame[rownames(GSE134108_frame) %in% miRNAs_to_plot, ]

expr_long_GSE134108 <- expr_subset_GSE134108 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)

metaplot_GSE134108 <- as.data.frame(meta_GSE134108_h)
metaplot_GSE134108$Sample <- rownames(metaplot_GSE134108)
rownames(metaplot_GSE134108) <- NULL

# Prepare predictions from results
predictions_GSE134108 <- data.frame(
  Sample = rownames(results_GSE134108),
  Prediction = results_GSE134108$max_score
)

head(predictions_GSE134108)


expr_long_GSE134108 <- expr_long_GSE134108 %>%
  left_join(metaplot_GSE134108, by = "Sample") %>%
  left_join(predictions_GSE134108, by = "Sample") %>%
  rename(actual_label = class)

metaplot_GSE134108 <- metaplot_GSE134108 %>%
  rename(actual_label = class)
predictions_GSE134108 <- predictions_GSE134108 %>%
  left_join(metaplot_GSE134108, by = "Sample")

prediction_outcomes_GSE134108 <- predictions_GSE134108 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))

head(prediction_outcomes_GSE134108)
table(prediction_outcomes_GSE134108$Correct)


# GSE140249_Validation ----------------------------------------------------

# Fetch the GEO data
GSE140249 <- getGEO("GSE140249", GSEMatrix = TRUE)
metadata_GSE140249 = pData(phenoData(GSE140249[[1]]))
metadata_GSE140249.subset = dplyr::select(metadata_GSE140249,c(2,34))
colnames(metadata_GSE140249.subset) <- c("Title", "Description")
frequency_table_GSE140249 <- table(metadata_GSE140249.subset$Description)

# organize
metadata.subset_cancer_GSE140249 = metadata_GSE140249.subset %>%
  dplyr::rename(state = "Description") %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))

str(GSE140249)
GSE140249_data <- GSE140249[[1]]
meta_GSE140249 <- pData(GSE140249_data)
annot_GSE140249 <- fData(GSE140249_data)
exp_GSE140249 <- exprs(GSE140249_data)
rownames(exp_GSE140249) <- annot_GSE140249$miRNA_ID_LIST

GSE140249_frame <- as.data.frame(exp_GSE140249)
GSE140249_frame <- na.omit(GSE140249_frame)

cleaned_GSE140249 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE140249_frame))
cleaned_GSE140249 <- gsub("-", "_", cleaned_GSE140249) 
cleaned_GSE140249 <- gsub(",", ".", cleaned_GSE140249) 
rownames(GSE140249_frame) <- cleaned_GSE140249
GSE140249mat <- as.matrix(GSE140249_frame)


meta_GSE140249_h <- read.csv("validation/GSE140249_metatest.csv", stringsAsFactors = TRUE,row.names = 1)
phenoData_GSE140249 <- new("AnnotatedDataFrame",data=meta_GSE140249_h)

GSE140249_expset <-ExpressionSet(GSE140249mat, phenoData = phenoData_GSE140249)

results_GSE140249 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                             Data = GSE140249_expset,
                                             tolerate_missed_genes = TRUE,
                                             weighted_votes = TRUE,
                                             verbose = TRUE)
knitr::kable(head(results_GSE140249))

confusion_GSE140249 <- caret::confusionMatrix(data = factor(results_GSE140249$max_score, 
                                                            levels = unique(train_object$data$Labels)),
                                              reference = factor(pData(GSE140249_expset)[,"class"], 
                                                                 levels = unique(train_object$data$Labels)),
                                              mode="everything")
print(confusion_GSE140249)
plot_confusion_matrix(confusion_GSE140249, "Confusion Matrix for GSE140249")

expr_subset_GSE140249 <- GSE140249_frame[rownames(GSE140249_frame) %in% miRNAs_to_plot, ]

expr_long_GSE140249 <- expr_subset_GSE140249 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)

metaplot_GSE140249 <- as.data.frame(meta_GSE140249_h)
metaplot_GSE140249$Sample <- rownames(metaplot_GSE140249)
rownames(metaplot_GSE140249) <- NULL

# Prepare predictions from results
predictions_GSE140249 <- data.frame(
  Sample = rownames(results_GSE140249),
  Prediction = results_GSE140249$max_score
)
head(predictions_GSE140249)

expr_long_GSE140249 <- expr_long_GSE140249 %>%
  left_join(metaplot_GSE140249, by = "Sample") %>%
  left_join(predictions_GSE140249, by = "Sample") %>%
  rename(actual_label = class)

metaplot_GSE140249 <- metaplot_GSE140249 %>%
  rename(actual_label = class)

# Merge the actual labels with predictions
predictions_GSE140249 <- predictions_GSE140249 %>%
  left_join(metaplot_GSE140249, by = "Sample")
prediction_outcomes_GSE140249 <- predictions_GSE140249 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))

head(prediction_outcomes_GSE140249)
table(prediction_outcomes_GSE140249$Correct)


# GSE150693_Validation ----------------------------------------------------

# Fetch the GEO data
GSE150693 <- getGEO("GSE150693", GSEMatrix = TRUE)

metadata_GSE150693 = pData(phenoData(GSE150693[[1]]))
metadata_GSE150693.subset = dplyr::select(metadata_GSE150693,c(1,37))
colnames(metadata_GSE150693.subset) <- c("Title", "Description")
frequency_table_GSE150693 <- table(metadata_GSE150693.subset$Description)

# organize
metadata.subset_cancer_GSE150693 = metadata_GSE150693.subset %>%
  dplyr::rename(state = "Description") %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))

str(GSE150693)
GSE150693_data <- GSE150693[[1]]
meta_GSE150693 <- pData(GSE150693_data)
annot_GSE150693 <- fData(GSE150693_data)
exp_GSE150693 <- exprs(GSE150693_data)
rownames(exp_GSE150693) <- annot_GSE150693$miRNA_ID_LIST

GSE150693_frame <- as.data.frame(exp_GSE150693)
GSE150693_frame <- na.omit(GSE150693_frame)

cleaned_GSE150693 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE150693_frame))
cleaned_GSE150693 <- gsub("-", "_", cleaned_GSE150693) 
cleaned_GSE150693 <- gsub(",", ".", cleaned_GSE150693) 
rownames(GSE150693_frame) <- cleaned_GSE150693
GSE150693mat <- as.matrix(GSE150693_frame)

meta_GSE150693_h <- read.csv("validation/GSE150693_metatest.csv", stringsAsFactors = TRUE,row.names = 1)
phenoData_GSE150693 <- new("AnnotatedDataFrame",data=meta_GSE150693_h)

GSE150693_expset <-ExpressionSet(GSE150693mat, phenoData = phenoData_GSE150693)

results_GSE150693 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                             Data = GSE150693_expset,
                                             tolerate_missed_genes = TRUE,
                                             weighted_votes = TRUE,
                                             verbose = TRUE)
knitr::kable(head(results_GSE150693))

confusion_GSE150693 <- caret::confusionMatrix(data = factor(results_GSE150693$max_score, 
                                                            levels = unique(train_object$data$Labels)),
                                              reference = factor(pData(GSE150693_expset)[,"class"], 
                                                                 levels = unique(train_object$data$Labels)),
                                              mode="everything")
print(confusion_GSE150693)
plot_confusion_matrix(confusion_GSE150693, "Confusion Matrix for GSE150693")

expr_subset_GSE150693 <- GSE150693_frame[rownames(GSE150693_frame) %in% miRNAs_to_plot, ]
expr_long_GSE150693 <- expr_subset_GSE150693 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)


metaplot_GSE150693 <- as.data.frame(meta_GSE150693_h)
metaplot_GSE150693$Sample <- rownames(metaplot_GSE150693)
rownames(metaplot_GSE150693) <- NULL

# Prepare predictions from results
predictions_GSE150693 <- data.frame(
  Sample = rownames(results_GSE150693),
  Prediction = results_GSE150693$max_score
)
head(predictions_GSE150693)

expr_long_GSE150693 <- expr_long_GSE150693 %>%
  left_join(metaplot_GSE150693, by = "Sample") %>%
  left_join(predictions_GSE150693, by = "Sample") %>%
  rename(actual_label = class)  # Ensure this matches your class label

metaplot_GSE150693 <- metaplot_GSE150693 %>%
  rename(actual_label = class)

# Merge the actual labels with predictions
predictions_GSE150693 <- predictions_GSE150693 %>%
  left_join(metaplot_GSE150693, by = "Sample")
prediction_outcomes_GSE150693 <- predictions_GSE150693 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))

head(prediction_outcomes_GSE150693)
table(prediction_outcomes_GSE150693$Correct)



# GSE85589_Validation -----------------------------------------------------

# Fetch the GEO data
GSE85589 <- getGEO("GSE85589", GSEMatrix = TRUE)
metadata_GSE85589 = pData(phenoData(GSE85589[[1]]))
metadata_GSE85589.subset = dplyr::select(metadata_GSE85589,c(2,8))
colnames(metadata_GSE85589.subset) <- c("Title", "Description")
frequency_table_GSE85589 <- table(metadata_GSE85589.subset$Description)

# organize
metadata.subset_cancer_GSE85589 = metadata_GSE85589.subset %>%
  dplyr::rename(state = "Description") %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))

str(GSE85589)
GSE85589_data <- GSE85589[[1]]
meta_GSE85589 <- pData(GSE85589_data)
annot_GSE85589 <- fData(GSE85589_data)
exp_GSE85589 <- exprs(GSE85589_data)
rownames(exp_GSE85589) <- annot_GSE85589$miRNA_ID
GSE85589_frame <- as.data.frame(exp_GSE85589)
GSE85589_frame <- na.omit(GSE85589_frame)

cleaned_GSE85589 <- gsub("hsa-miR-|hsa-let-", "", rownames(GSE85589_frame))
cleaned_GSE85589 <- gsub("-", "_", cleaned_GSE85589) 
cleaned_GSE85589 <- gsub(",", ".", cleaned_GSE85589) 
rownames(GSE85589_frame) <- cleaned_GSE85589
GSE85589mat <- as.matrix(GSE85589_frame)

meta_GSE85589_h <- read.csv("validation/metadata85589test.csv", stringsAsFactors = TRUE,row.names = 1)
phenoData_GSE85589 <- new("AnnotatedDataFrame",data=meta_GSE85589_h)

GSE85589_expset <-ExpressionSet(GSE85589mat, phenoData = phenoData_GSE85589)

results_GSE85589 <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                            Data = GSE85589_expset,
                                            tolerate_missed_genes = TRUE,
                                            weighted_votes = TRUE,
                                            verbose = TRUE)
knitr::kable(head(results_GSE85589))

confusion_GSE85589 <- caret::confusionMatrix(data = factor(results_GSE85589$max_score, 
                                                           levels = unique(train_object$data$Labels)),
                                             reference = factor(pData(GSE85589_expset)[,"class"], 
                                                                levels = unique(train_object$data$Labels)),
                                             mode="everything")
print(confusion_GSE85589)
plot_confusion_matrix(confusion_GSE85589, "Confusion Matrix for GSE85589")

expr_subset_GSE85589 <- GSE85589_frame[rownames(GSE85589_frame) %in% miRNAs_to_plot, ]

expr_long_GSE85589 <- expr_subset_GSE85589 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)

metaplot_GSE85589 <- as.data.frame(meta_GSE85589_h)
metaplot_GSE85589$Sample <- rownames(metaplot_GSE85589)
rownames(metaplot_GSE85589) <- NULL

# Prepare predictions from results
predictions_GSE85589 <- data.frame(
  Sample = rownames(results_GSE85589),
  Prediction = results_GSE85589$max_score
)

head(predictions_GSE85589)


expr_long_GSE85589 <- expr_long_GSE85589 %>%
  left_join(metaplot_GSE85589, by = "Sample") %>%
  left_join(predictions_GSE85589, by = "Sample") %>%
  rename(actual_label = class)

metaplot_GSE85589 <- metaplot_GSE85589 %>%
  rename(actual_label = class)

# Merge the actual labels with predictions
predictions_GSE85589 <- predictions_GSE85589 %>%
  left_join(metaplot_GSE85589, by = "Sample")

prediction_outcomes_GSE85589 <- predictions_GSE85589 %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))

head(prediction_outcomes_GSE85589)
table(prediction_outcomes_GSE85589$Correct)

# 5A_Validation_Datasets_Plot_Distribution -----------------------------------

metadesc_GSE120584 <- metadata.subset_cancer_GSE120584
metadesc_GSE120584$class <- meta_GSE120584_h$class

metadesc_GSE120584$description <- metadesc_GSE120584$state

# Create a named vector for replacements for shorter better visualization
replacement_map_GSE120584 <- c(
  "AD" = "Alzheimer's disease",
  "VaD" = "Vascular Dementia",
  "DLB" = "Dementia with Lewy Bodies",
  "MCI" = "Mild Cognitive Impairment",
  "NC" = "Negative Control"
)

# Add the description column and perform replacements
metadesc_GSE120584 <- metadesc_GSE120584 %>%
  mutate(description = recode(state, !!!replacement_map_GSE120584))

head(metadesc_GSE120584, 10)

metadesc_GSE117064 <- metadata_GSE117064.subset_filt
metadesc_GSE117064$class <- meta_GSE117064_h_filtered_df$class

metadesc_GSE117064$description <- metadesc_GSE117064$state

# Create a named vector for replacements
replacement_map_GSE117064 <- c(
  "CVD patient" = "cerebrovascular disorder (CVD)"
)

# Add the description column and perform replacements
metadesc_GSE117064 <- metadesc_GSE117064 %>%
  mutate(description = recode(Description, !!!replacement_map_GSE117064))

colnames(metadesc_GSE117064) <- colnames(metadesc_GSE120584)

head(metadesc_GSE117064, 10)


metadesc_GSE150693 <- metadata.subset_cancer_GSE150693
metadesc_GSE150693$class <- meta_GSE150693_h$class
metadesc_GSE150693$description <- metadesc_GSE150693$state

# Create a named vector for replacements
replacement_map_GSE150693 <- c(
  "MCI-NC" = "MCI-AD_Non_converted", "MCI-C" = "MCI-AD_Converted"
)

# Add the description column and perform replacements
metadesc_GSE150693 <- metadesc_GSE150693 %>%
  mutate(description = recode(state, !!!replacement_map_GSE150693))
head(metadesc_GSE150693, 10)


metadesc_GSE140249 <- metadata.subset_cancer_GSE140249
metadesc_GSE140249$class <- meta_GSE140249_h$class

metadesc_GSE140249$description <- metadesc_GSE140249$state

# Create a named vector for replacements
replacement_map_GSE140249 <- c(
  "AIH" = "Autoimmune_Hepatitis", "PBC" = "Primary_Biliary_Cholangitis", "OS" = "Overlap_Syndrome"
)

# Add the description column and perform replacements
metadesc_GSE140249 <- metadesc_GSE140249 %>%
  mutate(description = recode(state, !!!replacement_map_GSE140249))
head(metadesc_GSE140249, 10)


metadesc_GSE110317 <- metadata.subset_cancer_GSE110317
metadesc_GSE110317$class <- meta_GSE110317_h$class

metadesc_GSE110317$description <- metadesc_GSE110317$state

head(metadesc_GSE110317, 10)

metadesc_GSE85589 <- metadata.subset_cancer_GSE85589
metadesc_GSE85589$class <- meta_GSE85589_h$class

metadesc_GSE85589$description <- metadesc_GSE85589$state

# Create a named vector for replacements
replacement_map_GSE85589 <- c(
  "serum miRNA (pancreatic cancer)" = "Pancreatic cancer", "serum miRNA (intrahepatic cholangiocarcinoma)" = "Intrahepatic cholangiocarcinoma", "serum miRNA (stomach cancer)" = "Stomach cancer", "serum miRNA (colorectal cancer)" = "Colorectal cancer", "serum miRNA (GIST)" = "Gastrointestinal tract cancer", "serum miRNA (cholelithiasis)" = "Cholelithiasis", "serum miRNA (healthy control)" = "Healthy control", "serum miRNA (intrahepatic cholangiocarcinoma): validation" = "Intrahepatic cholangiocarcinoma")

# Add the description column and perform replacements
metadesc_GSE85589 <- metadesc_GSE85589 %>%
  mutate(description = recode(state, !!!replacement_map_GSE85589))

# View the first 10 rows of the updated data
head(metadesc_GSE85589, 10)

metadesc_GSE110651 <- metadata.subset_cancer_GSE110651
metadesc_GSE110651$class <- meta_GSE110651_h$class

metadesc_GSE110651$description <- metadesc_GSE110651$state

metadesc_GSE134108 <- metadata.subset_cancer_GSE134108
metadesc_GSE134108$class <- meta_GSE134108_h$class

metadesc_GSE134108$description <- metadesc_GSE134108$state

# Create a named vector for replacements
replacement_map_GSE134108 <- c(
  "Metastasis" = "Brain_cancer_metastasis", "no-Metastasis" = "Brain_cancer_no metastasis")

# Add the description column and perform replacements
metadesc_GSE134108 <- metadesc_GSE134108 %>%
  mutate(description = recode(state, !!!replacement_map_GSE134108))

# View the first 10 rows of the updated data
head(metadesc_GSE134108, 10)



# Combine all datasets into one data frame
datasets_val <- list(
  metadesc_GSE110317, metadesc_GSE85589, metadesc_GSE110651, 
  metadesc_GSE134108, metadesc_GSE120584, metadesc_GSE117064, 
  metadesc_GSE150693, metadesc_GSE140249
)

# Add a 'dataset' column to identify each dataset
names(datasets_val) <- c("GSE110317", "GSE85589", "GSE110651", 
                         "GSE134108", "GSE120584", "GSE117064", 
                         "GSE150693", "GSE140249")

# Combine datasets into a single data frame
combined_data_val <- do.call(rbind, lapply(names(datasets_val), function(name) {
  data <- datasets_val[[name]]
  data$dataset <- name
  return(data)
}))

# Summarize the counts for each state by class and dataset
library(dplyr)
summary_data_val <- combined_data_val %>%
  group_by(dataset, description, class) %>%
  summarise(count = n(), .groups = "drop")

# Load ggplot2 for visualization
library(ggplot2)

# Create the geom_tile plot
ggplot(summary_data_val, aes(x = dataset, y = description, fill = class)) +
  geom_tile(color = "white") +
  geom_text(aes(label = count), color = "white", size = 6, fontface = "bold") +  # Bold text for counts
  scale_fill_manual(values = c("cancer" = "red", "non_cancer" = "blue"), name = "Class") +
  labs(title = "Distribution of Samples in Validation datasets",
       x = "Dataset", 
       y = "State") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black", size = 12),  # Bold x-axis labels
    axis.text.y = element_text(face = "bold", color = "black", size = 12),                        # Bold y-axis labels
    axis.title = element_text(face = "bold", size = 16),                         # Bold axis titles
    plot.title = element_text(face = "bold", hjust = 0.5),            # Bold and centered title
    legend.title = element_text(face = "bold"),                       # Bold legend title
    legend.text = element_text(face = "bold", size = 16)                         # Bold legend text
  )

ggsave("../results/Figures/5A_validationdistribution_reduced.png", dpi = 600,  height = 8, width = 10)


# 5B_Confusion Matrices_Validation -------------------------------------------

#Confusion_Matrices_Validation

# Sample confusion matrices in list form for demonstration
confusion_matrices_validation <- list(
  GSE110317 = confusion_GSE110317,
  GSE110651 = confusion_GSE110651,
  GSE117064 = confusion_GSE117064,
  GSE120584 = confusion_GSE120584,
  GSE134108 = confusion_GSE134108,
  GSE140249 = confusion_GSE140249,
  GSE150693 = confusion_GSE150693,
  GSE85589 = confusion_GSE85589)

# Convert each matrix to a data frame and add a `Dataset` column
all_cm_validation <- do.call(rbind, lapply(names(confusion_matrices_validation), function(name) {
  cm <- as.data.frame(confusion_matrices_validation[[name]]$table)  # assumes confusion_matrix$table format
  cm$Dataset <- name  # add a column to indicate the dataset
  cm  # return modified data frame
}))


ggplot(all_cm_validation, aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), color = "black", fontface = "bold", size = 7) +  # Bold the text in tiles
  facet_wrap(~ Dataset, ncol = 2) +
  scale_fill_gradient(low = "lavender", high = "lightblue") +
  theme_minimal() +
  labs(title = "Confusion Matrices - Test Datasets",
       x = "Predicted Label", y = "True Label") +
  theme(
    plot.title = element_text(face = "bold", size = 8, hjust = 0.5),  # Bold and center the title
    axis.title.x = element_text(face = "bold", size = 16, color = "black"),  # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 16, color = "black"),  # Bold y-axis title
    axis.text.x = element_text(face = "bold", size = 16, color = "black"),   # Bold x-axis text
    axis.text.y = element_text(face = "bold", size = 16, color = "black"),   # Bold y-axis text
    strip.text = element_text(face = "bold", size = 20, color = "black"),    # Bold facet titles
    legend.text = element_text(face = "bold", size = 16),  # Bold and increase legend text size
    legend.title = element_text(face = "bold", size = 16)  # Bold and increase legend title size
  )

ggsave("../results/Figures/5B_Confusion_Matrix_withValidation.png", dpi = 600, height = 8, width = 8)

# 5D_Incorrect_Validation -------------------------------------------------

validation_datasets <- c("GSE120584", "GSE117064", "GSE150693", "GSE140249", "GSE110651", "GSE110317","GSE85589")

incorrect_samples_validation <- list()

# Loop over each validation dataset
for (dataset in validation_datasets) {
  
  # Load the prediction outcomes for the dataset
  prediction_outcomes <- get(paste0("prediction_outcomes_", dataset))
  
  # Filter out incorrectly predicted samples
  incorrect_samples <- prediction_outcomes[prediction_outcomes$Correct == "Incorrect", ]
  
  # Load the corresponding metadata for the dataset
  meta_data <- get(paste0("metadesc_", dataset))
  
  # Ensure the Sample column is a character type and remove any leading/trailing spaces
  incorrect_samples$Sample <- trimws(as.character(incorrect_samples$Sample))
  meta_data$class <- trimws(as.character(meta_data$class))  # Ensure the "class" column is trimmed
  
  # Rename 'class' column to 'Sample' in meta_data to match the 'prediction_outcomes' data frame
  colnames(meta_data)[colnames(meta_data) == "class"] <- "Sample"
  
  # Convert row names of meta_data to the 'MetaSample' column to avoid duplicate 'Sample' column
  meta_data <- rownames_to_column(meta_data, var = "MetaSample")
  
  # Now, rename 'MetaSample' to 'Sample' (if you still want to retain this name)
  meta_data$Sample <- meta_data$MetaSample
  meta_data <- meta_data %>% dplyr::select(-MetaSample)  # Explicitly use dplyr::select to avoid conflicts
  
  # Merge the incorrect samples with the metadata
  incorrect_samples_meta <- merge(incorrect_samples, meta_data, by = "Sample")
  
  # Add a column for the dataset name
  incorrect_samples_meta$dataset <- dataset
  
  # Append the result to the list
  incorrect_samples_validation[[dataset]] <- incorrect_samples_meta
}

# Combine all the validation datasets into one data frame
incorrect_samples_validation_combined <- do.call(rbind, incorrect_samples_validation)

# Process the combined validation data: Calculate counts and color
incorrect_samples_validation_combined <- incorrect_samples_validation_combined %>%
  group_by(dataset, description) %>%
  mutate(count = n()) %>%  # Calculate the count for each combination
  ungroup() %>%
  mutate(
    class_label = ifelse(description == "cancer", "Cancer", "Non-Cancer"),  # Class labels for color
    class_color = ifelse(description == "cancer", "red", "blue"),  # Color for class (Cancer - red, Non-Cancer - green)
    fill_color = case_when(
      actual_label == "cancer" & Prediction == "non_cancer" ~ "red",  # Highlight with red for incorrect cancer prediction
      count == 0 ~ "lightgray",  # Mild gray for zero counts
      TRUE ~ class_color  # Default to class_color for other cases
    )
  )

# Create the tile plot for validation datasets
ggplot(incorrect_samples_validation_combined, aes(x = dataset, y = description, fill = fill_color)) +
  geom_tile() +  # Use tiles for better readability
  geom_text(aes(label = count), size = 6, color = "white") +  # Display count inside each tile
  scale_fill_identity() +  # Use the fill_color directly
  labs(
    title = "Distributions of Incorrect Samples in Validation Datasets",
    x = "Dataset",
    y = "State",
    fill = "Sample Class"
  ) +
  theme_minimal() +
  theme(
    # Ensure the Y-axis labels are black
    axis.text.y = element_text(
      size = 15, 
      face = "bold", 
      color = "black"  # Set all Y-axis labels to black
    ),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold", color = "black"),
    strip.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(face = "bold", size = 16),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  )
# Save the plot for validation datasets
ggsave("../results/Figures/5D_incorrectplots_validation_benign_marked.png", height = 8, width = 14)

# 5F_Test_IncorrectPlots -----------------------------------------------------

# Initialize the list to store the dataset results
dataset_list <- list()

# Define the datasets
datasets <- c("GSE106817", "GSE112264", "GSE113740", 
              "GSE122497", "GSE137140", "GSE139031", "GSE164174", "GSE73002")

# Loop through each dataset and prepare the data
for (dataset in datasets) {
  # Assuming metafull and prediction outcomes are available for each dataset
  metafull <- get(paste0("metafull_", dataset))
  prediction_outcomes <- get(paste0("prediction_outcomes_", dataset))
  
  # Convert rownames to column (assuming 'Sample' is the rowname)
  metafull <- rownames_to_column(metafull, var = "Sample")
  
  # Rename 'description' column in metafull to avoid conflict with prediction_outcomes
  metafull <- metafull %>%
    rename(description_metafull = description)
  
  # Merge data: match samples and add description for those samples
  merged_data <- merge(prediction_outcomes, metafull, by = "Sample", all.x = TRUE)
  
  # Check if the 'description' column is now in merged_data
  print(paste("Columns in merged data for dataset", dataset, ":", colnames(merged_data)))
  
  # Filter out only the incorrect predictions using dplyr::select to avoid conflicts
  incorrect_data <- merged_data %>%
    filter(Correct != "Correct") %>%
    dplyr::select(Sample, description_metafull, Correct) %>%  # Use renamed description column
    mutate(Dataset = dataset)  # Add the dataset column
  
  # Append to the list of dataset data
  dataset_list[[dataset]] <- incorrect_data
}

# Combine all datasets into one data frame
combined_data_inc <- bind_rows(dataset_list)

# Plotting the data using ggplot2
ggplot(combined_data_inc, aes(x = Dataset, y = Sample, fill = Correct)) +
  geom_tile() +
  scale_fill_manual(values = c("Incorrect" = "red", "Correct" = "green")) + 
  labs(title = "Incorrect Predictions Across Datasets", 
       x = "Dataset", 
       y = "Sample ID", 
       fill = "Prediction Status") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6),  # Reduce font size for Sample labels
        axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x labels for clarity



# Initialize a list to store the incorrect samples for each dataset
incorrect_samples_all <- list()

# Loop over each dataset
for (dataset in datasets) {
  
  # Load the prediction outcomes for the dataset
  prediction_outcomes <- get(paste0("prediction_outcomes_", dataset))
  
  # Filter out incorrectly predicted samples
  incorrect_samples <- prediction_outcomes[prediction_outcomes$Correct == "Incorrect", ]
  
  # Load the corresponding metadata for the dataset
  meta_data <- get(paste0("metafull_", dataset))
  
  # Ensure the Sample column is a character type and remove any leading/trailing spaces
  incorrect_samples$Sample <- trimws(as.character(incorrect_samples$Sample))
  meta_data$class <- trimws(as.character(meta_data$class))  # Ensure the "class" column is trimmed
  
  # Rename 'class' column to 'Sample' in meta_data to match the 'prediction_outcomes' data frame
  colnames(meta_data)[colnames(meta_data) == "class"] <- "Sample"
  
  # Convert row names of meta_data to the 'MetaSample' column to avoid duplicate 'Sample' column
  meta_data <- rownames_to_column(meta_data, var = "MetaSample")
  
  # Now, rename 'MetaSample' to 'Sample' (if you still want to retain this name)
  meta_data$Sample <- meta_data$MetaSample
  
  # Drop the temporary MetaSample column using base R
  meta_data <- meta_data[, !names(meta_data) %in% "MetaSample"]
  
  # Merge the incorrect samples with the metadata
  incorrect_samples_meta <- merge(incorrect_samples, meta_data, by = "Sample")
  
  # Add a column for the dataset name
  incorrect_samples_meta$dataset <- dataset
  
  # Append the result to the list
  incorrect_samples_all[[dataset]] <- incorrect_samples_meta
}

# Combine all the datasets into one data frame
incorrect_samples_all_combined <- do.call(rbind, incorrect_samples_all)

incorrect_samples_all_combined <- incorrect_samples_all_combined %>%
  group_by(dataset, description.x) %>%  # Use description.x instead of description
  mutate(count = n()) %>%  # Calculate the count for each combination
  ungroup() %>%
  mutate(
    class_label = ifelse(description.x == "cancer", "Cancer", "Non-Cancer"),  # Class labels for color
    class_color = ifelse(description.x == "cancer", "red", "blue"),  # Color for class (Cancer - red, Non-Cancer - blue)
    fill_color = case_when(
      actual_label == "cancer" & Prediction == "non_cancer" ~ "red",  # Highlight with red for incorrect cancer prediction
      count == 0 ~ "lightgray",  # Mild gray for zero counts
      TRUE ~ class_color  # Default to class_color for other cases
    )
  )

# Process the combined data: Calculate counts and color
incorrect_samples_all_combined <- incorrect_samples_all_combined %>%
  group_by(dataset, description.x) %>%
  mutate(count = n()) %>%  # Calculate the count for each combination
  ungroup() %>%
  mutate(
    class_label = ifelse(description.x == "cancer", "Cancer", "Non-Cancer"),  # Class labels for color
    class_color = ifelse(description.x == "cancer", "red", "blue"),  # Color for class (Cancer - red, Non-Cancer - green)
    fill_color = case_when(
      actual_label == "cancer" & Prediction == "non_cancer" ~ "red",  # Highlight with red for incorrect cancer prediction
      count == 0 ~ "lightgray",  # Mild gray for zero counts
      TRUE ~ class_color  # Default to class_color for other cases
    )
  )

# Create the tile plot with proper Y-axis label coloring and zero tile colors
ggplot(incorrect_samples_all_combined, aes(x = dataset, y = description.x, fill = fill_color)) +
  geom_tile() +  # Use tiles for better readability
  geom_text(aes(label = count), size = 5, color = "white") +  # Display count inside each tile
  scale_fill_identity() +  # Use the fill_color directly
  labs(
    title = "Distributions of Incorrect Samples in Test Datasets",
    x = "Dataset",
    y = "State",
    fill = "Sample Class"
  ) +
  theme_minimal() +
  theme(
    # Ensure the Y-axis labels are black
    axis.text.y = element_text(
      size = 12, 
      face = "bold", 
      color = "black"  # Set all Y-axis labels to black
    ),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", size = 16),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave("../results/Figures/5F_incorrectpplots_benign_marked.png", dpi = 600, height = 12, width = 18)

# 5H_Incorrect_Train ------------------------------------------------------

meta_train_inc <- meta_train
meta_train_inc$description <- metadata.subset_GSE211692$state

meta_train_inc[1:8,1:2]

metafull <- rownames_to_column(meta_train_inc, var = "Sample")
meta_train_inc <- rownames_to_column(meta_train_inc, var = "Sample")
meta_train_inc$Sample <- trimws(meta_train_inc$Sample)

metaplot_train <- as.data.frame(meta_train)
metaplot_train$Sample <- rownames(meta_train)
rownames(metaplot_train) <- NULL

# Prepare predictions from results
predictions_train <- data.frame(
  Sample = rownames(results_train),  # Use row names as Sample names
  Prediction = results_train$max_score
)
head(predictions_train)

metaplot_train <- metaplot_train %>%
  rename(actual_label = class)

# Merge the actual labels with predictions
predictions_train <- predictions_train %>%
  left_join(metaplot_train, by = "Sample")

prediction_outcomes_train <- predictions_train %>%
  mutate(Correct = ifelse(actual_label == Prediction, "Correct", "Incorrect"))
head(prediction_outcomes_train)
table(prediction_outcomes_train$Correct)

prediction_outcomes_train$Sample <- trimws(prediction_outcomes_train$Sample)

# Merge prediction outcomes with metadata using 'Sample' as the key
merged_data <- merge(prediction_outcomes_train, meta_train_inc, by = "Sample", all.x = TRUE)

head(merged_data)


incorrect_data <- merged_data[merged_data$Correct != "Correct", c("Sample", "description", "Correct")]


# Calculate counts for each description
description_distribution <- merged_data %>%
  filter(Correct != "Correct") %>%
  group_by(description) %>%
  summarize(count = n()) %>%
  ungroup()


# Calculate percentage for labels
description_distribution <- description_distribution %>%
  mutate(percentage = (count / sum(count)) * 100)

description_distribution_split <- description_distribution %>%
  mutate(description = ifelse(description == "extraparenchymal brain tumor and benign brain", 
                              "extraparenchymal brain tumor\nand benign brain", 
                              description))

ggplot(description_distribution_split, aes(x = count, y = reorder(description, count))) +
  geom_bar(stat = "identity", fill = "blue", color = "black") +  # Uniform bar color
  geom_text(aes(label = count), 
            hjust = -0.2, size = 3.5, fontface = "bold", color = "black") +  # Dark text and number-only labels
  labs(title = "Incorrect predictions_Training_data",
       x = "Count",
       y = "Description") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # Bold centered title
    axis.title.x = element_text(size = 14, face = "bold"),  # Bold x-axis label
    axis.title.y = element_text(size = 12, face = "bold"),  # Bold y-axis label
    axis.text.x = element_text(size = 12, face = "bold"),  # Bold x-axis tick labels
    axis.text.y = element_text(size = 12, face = "bold")   # Bold y-axis tick labels (category labels)
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 350))  # Reduce the x-axis length


ggsave("../results/Figures/5H_incorrectdistribution_train.png", dpi = 600, height = 8, width = 12)


# Validation_merged -------------------------------------------------------

# Validation Datasets
validation_datasets <- list(
  GSE110317 = GSE110317_expset,
  GSE110651 = GSE110651_expset,
  GSE117064 = GSE117064_expset_filtered,
  GSE120584 = GSE120584_expset,
  GSE134108 = GSE134108_expset,
  GSE140249 = GSE140249_expset,
  GSE150693 = GSE150693_expset,
  GSE85589  = GSE85589_expset
)

#Find common features (genes) across all validation datasets
common_features_val <- Reduce(intersect, lapply(validation_datasets, function(x) featureNames(x)))
validation_datasets_subset <- lapply(validation_datasets, function(x) x[common_features_val, ])

# Merge expression matrices
expression_matrices_val <- lapply(validation_datasets_subset, exprs)
merged_exprs_val <- do.call(cbind, expression_matrices_val)

# Merge sample metadata and include dataset name
metadata_list_val <- lapply(names(validation_datasets_subset), function(dataset_name) {
  meta <- pData(validation_datasets_subset[[dataset_name]])
  meta$Dataset <- dataset_name  # Add dataset name as a new column
  return(meta)
})

merged_metadata_val <- do.call(rbind, metadata_list_val)
merged_meta_val <- merged_metadata_val
merged_meta_val$Dataset <- NULL 
merged_validation_mat <- as.matrix(merged_exprs_val)

phenoData_mergedval <- new("AnnotatedDataFrame", data = merged_meta_val)
merged_expressionset_val <- ExpressionSet(merged_validation_mat, phenoData = phenoData_mergedval)

results_validation_merged <- predict_one_vs_rest_TSP(
  classifier = classifier_train,
  Data = merged_expressionset_val,
  tolerate_missed_genes = TRUE,
  weighted_votes = TRUE,
  verbose = TRUE
)
knitr::kable(head(results_validation_merged))

confusion_validation_merged <- caret::confusionMatrix(
  data = factor(results_validation_merged$max_score, 
                levels = unique(train_object$data$Labels)),
  reference = factor(pData(merged_expressionset_val)[,"class"], 
                     levels = unique(train_object$data$Labels)),
  mode = "everything"
)
print(confusion_validation_merged)
plot_confusion_matrix(confusion_validation_merged, "Confusion Matrix for all Validation Datasets")

ggsave("../results/Figures/confusion_validation_mergeddatasets.png")


descrip_val <- combined_data_val[, c("class", "description", "dataset")]


unique_descriptions_val <- unique(descrip_val$description)
benign_keyword <- "Benign"
benign_color <- "#A3A500"
num_non_benign_val <- sum(!sapply(unique_descriptions_val, function(x) grepl(benign_keyword, x, ignore.case = TRUE)))

custom_palette <- brewer.pal(12, "Set3")  # Choose a palette with no red/blue dominance
custom_palette <- custom_palette[!custom_palette %in% c("#FB8072", "#80B1D3")]  # Exclude specific red/blue shades

if (num_non_benign_val > length(custom_palette)) {
  custom_palette <- colorRampPalette(custom_palette)(num_non_benign_val)
}

# Initialize colors vector
colors_val <- character(length(unique_descriptions_val))

# Assign colors
color_index_val <- 1
for (i in seq_along(unique_descriptions_val)) {
  if (grepl(benign_keyword, unique_descriptions_val[i], ignore.case = TRUE)) {
    colors_val[i] <- benign_color
  } else {
    colors_val[i] <- custom_palette[color_index_val]
    color_index_val <- color_index_val + 1
  }
}

# Create a named vector for clarity
names(colors_val) <- unique_descriptions_val

# Display assigned colors
colors_val


descrip_val <- descrip_val %>%
  mutate(disease_types = case_when(
    grepl("cancer|Intrahepatic cholangiocarcinoma", description, ignore.case = TRUE) ~ "cancer",
    grepl("benign|Benign", description, ignore.case = TRUE) ~ "benign",
    grepl("Alzheimer|Dementia|MCI|Cognitive|Vascular|Mild Cognitive Impairment", description, ignore.case = TRUE) ~ "Neurological diseases",
    grepl("Autoimmune|Hepatitis|Primary_Biliary|Cholangitis|Overlap", description, ignore.case = TRUE) ~ "autoimmune diseases",
    grepl("Healthy|Healthy control|Non-CVD control|Negative Control", description, ignore.case = TRUE) ~ "Healthy control",
    TRUE ~ "Other"
  ))

head(descrip_val, 10)


color_map <- c(
  "cancer" = "#FF5733",  # Red
  "benign" = "#A3A500",  # Green
  "Neurological diseases" = "#3498DB",  # Blue
  "autoimmune diseases" = "#F39C12",  # Orange
  "Healthy control" = "#2ECC71",  # Light Green
  "Other" = "#28B463"  # Gray
)

ref_colors <- c("cancer" = "red", "non_cancer" = "blue")

png('../results/Figures/validation_merged.png', res = 600, units = "in", width = 16, height = 10, bg = "white")

# Plot the heatmap with custom reference label colors
plot_binary_TSP(Data = merged_expressionset_val,               # Your data object
                classifier = classifier_train,     # Your classifier
                ref = "class",
                binary_col = c("salmon", "lightgreen", "gray"), # Binary heatmap colors
                prediction = results_validation_merged,        # Your prediction data
                platform = descrip_val$description,  # Platform/study info
                platform_col = colors_val,     # Platform colors
                show_rule_name = TRUE,        # Show rule names
                legend = FALSE,               # Hide legend
                anno_height = 0.04,          # Annotation height
                score_height = 0.075,         # Score height
                title = "cancer_classifier_Validation",  # Title of the plot
                ref_col = ref_colors,
                pred_col = ref_colors,
                margin = c(0,6,0,6))
dev.off()


#platform_legend_val <- Legend(
  #title = "Platform/Stud",
  #at = names(colors_val),
  #legend_gp = gpar(fill = colors_val),   # Colors for the legend
  #ncol = 1 ,                                  # Number of rows in the legend
  #title_gp = gpar(fontsize = 14, fontface = "bold", col = "black"),  # Bold, black title
  #labels_gp = gpar(fontsize = 12, fontface = "bold", col = "black")  # Bold, black labels
#)


# Save the platform legend as a separate plot
#png('../results/Figures/platform_legend_val.png', res = 600, units = "in", width = 14, height = 6, bg = "white")
#draw(platform_legend_val)
#dev.off()
