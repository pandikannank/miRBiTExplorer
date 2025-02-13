#This code will reproduce the pre-processing of Training Data and Training the classifier
#Fig-2A,2C,2D,2E,2F,2G,2H,3A,3B,3D,3E
########################Preprocessing###########################################


#Training_data_GSE211692
training_data = read.delim('train/GSE211692_processed_data_matrix.txt', header = T)

#Download the clinical metadata from GEO
GSE211692= getGEO(GEO = "GSE211692",GSEMatrix = TRUE)
metadata_GSE211692 = pData(phenoData(GSE211692[[1]]))

#Harmonize column names into GEO_sample_IDs_for consistency and uniqueness
title_to_geo_map <- metadata_GSE211692$geo_accession
names(title_to_geo_map) <- metadata_GSE211692$title

# Replace column names in training_data
colnames(training_data) <- sapply(colnames(training_data), function(col) {
  if (col %in% names(title_to_geo_map)) {
    return(title_to_geo_map[col])
  } else {
    return(col)
  }
})

#select the necessary columns
metadata.subset_GSE211692 = dplyr::select(metadata_GSE211692,c(2,10))
frequency_table_GSE211692 <- table(metadata.subset_GSE211692$characteristics_ch1)


# 2A_Disease state Visualization ------------------------------------------

metadisttrain <- metadata.subset_GSE211692
metadisttrain$class <- meta_train$class
metadisttrain$title <- rownames(meta_train)

# Prepare data for plotting
state_data <- metadisttrain %>%
  count(state, class) %>%
  arrange(state, class)  # Sort by state first, then by class to make sure cancer comes first

# Convert state to a factor to control the order on y-axis
state_data$state <- factor(state_data$state, levels = unique(state_data$state))

ggplot(state_data, aes(x = n, y = state, color = class, group = state)) +
  geom_segment(aes(x = 0, xend = n, y = state, yend = state), color = "grey80", size = 0.5) +  # Background lines for each state
  geom_point(aes(color = class), size = 6, position = position_dodge(width = 0.5)) +  # Offset points for clarity
  geom_text(aes(label = n), hjust = -0.5, size = 4, position = position_dodge(width = 0.5)) +  # Offset labels for clarity
  scale_color_manual(values = c("cancer" = "red", "non_cancer" = "blue")) +  # Custom colors for classes
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16)
  ) +
  labs(
    title = "Training data Distribution - GSE211692",
    x = "Count",
    y = "State",
    color = "Class"
  ) +
  coord_cartesian(clip = "off")  # Ensure labels outside the plot boundary are visible

ggsave("Figures/2ASample_distribution_plot.png",height = 8, width = 10)

# organize
metadata.subset_GSE211692 = metadata_GSE211692 %>%
  dplyr::select(2,10) %>%
  dplyr::rename(state = characteristics_ch1 ) %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))

# reshaping the data for further operations
dat.long_GSE211692 = training_data %>%
  gather(key= "sample", value = "expression", -ID_REF)


# join 
dat.long_GSE211692 = dat.long_GSE211692 %>%
  left_join(.,metadata.subset_GSE211692,by = c("sample" = "geo_accession"))

cancervsnoncancer = dat.long_GSE211692[,-4]

# Reshape data to wide format
cvsnc <- cancervsnoncancer %>% spread(key = sample, value = expression)

# Set row names and remove ID_REF column
row.names(cvsnc) <- cvsnc$ID_REF
cvsnc <- cvsnc[,-1]

# Filter out low-expression miRNAs
min_expr <- 5
keep <- rowMeans(cvsnc) > min_expr
cvsnc_filtered <- cvsnc[keep, ]

# Log2 transformation, if required
qx <- as.numeric(quantile(cvsnc_filtered, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) { 
  cvsnc_filtered[which(cvsnc_filtered <= 0)] <- NaN
  cvsnc_filtered <- log2(cvsnc_filtered) 
}

# Normalize the filtered data
Nx_filtered <- normalizeBetweenArrays(cvsnc_filtered)

# Check for NAs presence
sum(is.na(Nx_filtered))

expression_quantiles <- apply(Nx_filtered, 1, quantile, probs = 0.50, na.rm = TRUE)
highly_expressed_train <- Nx_filtered[expression_quantiles > quantile(expression_quantiles, 0.50), ]
print(dim(highly_expressed_train))

#miRNA_names_shortening for better simplicity and better visualization in plots
cleanedDatah <- gsub("hsa-miR-|hsa-let-", "", rownames(highly_expressed_train))
cleanedDatah <- gsub("-", "_", cleanedDatah) 
cleanedDatah <- gsub(",", ".", cleanedDatah) 
rownames(highly_expressed_train) <- cleanedDatah
highframe_train <- as.data.frame(highly_expressed_train)


# 2C_t-SNE_Visualization --------------------------------------------------

umap_result <- umap(t(highframe_train))# Plot UMAP with specific colors for cancer and non-cancer
umap_coords <- as.data.frame(umap_result$layout)
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$State <- meta_train$class[match(rownames(umap_coords), metadata.subset_GSE211692$geo_accession)]
umap_coords$Category <- ifelse(umap_coords$State == "cancer", "cancer", "non_cancer")
unique(umap_coords$State)

tsne_result <- Rtsne(t(highly_expressed_train), dims = 2, pca = TRUE, check_duplicates = FALSE)

# Extract t-SNE coordinates
tsne_coords <- as.data.frame(tsne_result$Y)
tsne_coords$Category <- umap_coords$Category

# Plot t-SNE
ggplot(tsne_coords, aes(x = V1, y = V2, color = Category)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Plot of Cancer vs. Non-Cancer",
       x = "t-SNE 1",
       y = "t-SNE 2") +
  scale_color_manual(values = c("cancer" = "red", "non_cancer" = "blue")) +
  theme(legend.position = "right")

ggsave("Figures/2Ct-sne_train.png")


###################Training Single Sample Classifier############################

meta_train <- read.csv("train/meta_train_full.csv", row.names = 1, stringsAsFactors = TRUE)
frequency_table_train <- table(meta_train$class)

phenoData_train<- new("AnnotatedDataFrame",data=meta_train)

#Convert to ExpressionSet format
train_high <- ExpressionSet(highly_expressed_train, phenoData = phenoData_train)

#Training
train_object <- ReadData(Data = train_high, 
                         Labels = "class", 
                         Platform = NULL, 
                         verbose = TRUE)
#optional
filtered_genes_train <- filter_genes_TSP(data_object = train_object,
                                         filter = "one_vs_rest",
                                         platform_wise = FALSE,
                                         UpDown = TRUE,
                                         verbose = TRUE)
filtered_genes_train

#Classifier_Training
classifier_train <- train_one_vs_rest_TSP(data_object = train_object,
                                          filtered_genes = filtered_genes_train,
                                          k_range = 2:5,
                                          include_pivot = TRUE,
                                          one_vs_one_scores = TRUE,
                                          platform_wise_scores = FALSE,
                                          seed = 1234,
                                          verbose = TRUE)
#results of the classifier
results_train <- predict_one_vs_rest_TSP(classifier = classifier_train,
                                         Data = train_object,
                                         tolerate_missed_genes = TRUE,
                                         weighted_votes = TRUE,
                                         verbose = TRUE)


# 2D_Rules_Visualization --------------------------------------------------

# Create the cancer rules data frame by directly accessing TSPs as a matrix
cancer_rules_train <- data.frame(
  Rule = rownames(classifier_train$classifiers$cancer$TSPs),
  Gene1 = classifier_train$classifiers$cancer$TSPs[, "gene1"],
  Gene2 = classifier_train$classifiers$cancer$TSPs[, "gene2"],
  Score = classifier_train$classifiers$cancer$score,
  TieVote = classifier_train$classifiers$cancer$tieVote
)

# View the created data frame
print(cancer_rules_train)

#write.csv(cancer_rules, "cancer_classifier_TSP_rules4.csv", row.names = FALSE)


# Create the cancer rules data frame by directly accessing TSPs as a matrix
no_rules_train <- data.frame(
  Rule = rownames(classifier_train$classifiers$no$TSPs),
  Gene1 = classifier_train$classifiers$no$TSPs[, "gene1"],
  Gene2 = classifier_train$classifiers$no$TSPs[, "gene2"],
  Score = classifier_train$classifiers$no$score,
  TieVote = classifier_train$classifiers$no$tieVote
)

# View the created data frame
print(no_rules_train)

#Rules for cancer
classifier_train$classifiers$cancer$TSPs

#Rules for non_cancer
classifier_train$classifiers$cancer$TSPs

# Extract rules for 'Cancer' class
cancer_rules <- classifier_train$classifiers$cancer$TSPs
cancer_scores <- classifier_train$classifiers$cancer$score

# Create a data frame for the 'cancer' class rules
cancer_df <- data.frame(
  miRNA_pair = apply(cancer_rules, 1, paste, collapse = ", "), 
  score = cancer_scores
)

# Create a data frame for the 'non_cancer' class rules
non_cancer_rules <- classifier_train$classifiers$non_cancer$TSPs
non_cancer_scores <- classifier_train$classifiers$non_cancer$score

# Create a data frame for the 'non_cancer' class
non_cancer_df <- data.frame(
  miRNA_pair = apply(non_cancer_rules, 1, paste, collapse = ", "), 
  score = non_cancer_scores
)

# Add a column to indicate the class
cancer_df$class <- "cancer"
non_cancer_df$class <- "non_cancer"

# Combine both data frames into one for plotting
combined_df <- rbind(cancer_df, non_cancer_df)

# Reorder the miRNA_pair within each class based on score
combined_df$miRNA_pair <- factor(combined_df$miRNA_pair, levels = unique(combined_df$miRNA_pair[order(combined_df$score)]))
combined_df$miRNA_pair <- factor(combined_df$miRNA_pair, levels = combined_df$miRNA_pair)

ggplot(combined_df, aes(x = miRNA_pair, y = score, fill = class)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # Adjust width for better spacing
  coord_flip() +  # Flip the coordinates for horizontal bars
  geom_text(aes(label = round(score, 3)), nudge_x = 0.2, hjust = -0.25, vjust =1, color = "black", fontface = "bold", size = 5) +  # Add space to the right of the bars for score labels
  scale_fill_manual(values = c("cancer" = "red", "non_cancer" = "blue")) +  # Color the bars differently for each class
  scale_y_continuous(limits = c(0, 1), breaks = seq(0.1, 1, 0.1)) +  # Set y-axis scale from 0.1 to 1 with breaks at 0.1 intervals
  theme_minimal(base_size = 12) +  # Base size for fonts
  labs(
    title = "miRNAs Pair Scores by Class", 
    x = "miRNAs Pair", 
    y = "Score"
  ) +
  theme(
    axis.text.y = element_text(size = 16, face = "bold", color = "black"),  # Bold and black y-axis labels
    axis.text.x = element_text(size = 16, face = "bold", color = "black"),  # Bold and black x-axis labels
    axis.title = element_text(face = "bold", color = "black", size = 14),  # Bold and black axis titles
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),  # Bold and centered plot title
    panel.grid.major = element_line(color = "gray", size = 0.2),  # Lighter gridlines for cleaner look
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add square border around the plot # Remove the box outline
    plot.margin = margin(10, 10, 10, 40)  # Add more space on the left and right side of the plot
  )

#ggsave("Figures/2D_mirna_pairs.png", height = 12, width = 12)


#Individual_Rules_Performance
actual_labels <- meta_train$class
top_rules_train <- cancer_rules_train %>%
  arrange(desc(Score)) %>%
  slice(1:5)
miRNAs_to_plot <- unique(c(top_rules_train$Gene1, top_rules_train$Gene2))
expr_subset_train <- highframe_train[rownames(highframe_train) %in% miRNAs_to_plot, ]

# Store metrics in a new data frame
rule_results <- data.frame(
  Rule = character(),
  Accuracy = numeric(),
  Sensitivity = numeric(),
  Specificity = numeric(),
  Precision = numeric(),
  F1 = numeric(),
  Recall = numeric(),
  stringsAsFactors = FALSE
)

# Iterate through each rule to calculate performance metrics
for (i in 1:nrow(cancer_rules_train)) {
  # Get miRNAs for the current rule
  gene1 <- cancer_rules_train$Gene1[i]
  gene2 <- cancer_rules_train$Gene2[i]
  rule_name <- cancer_rules_train$Rule[i]
  
  # Extract expression data for the two miRNAs
  expr_gene1 <- expr_subset_train[gene1, ]
  expr_gene2 <- expr_subset_train[gene2, ]
  
  # Apply the rule (e.g., if miRNA 1 > miRNA 2, predict "cancer", else predict "non_cancer")
  predictions <- ifelse(expr_gene1 > expr_gene2, "cancer", "non_cancer")
  
  # Create confusion matrix
  cm <- confusionMatrix(factor(predictions, levels = c("cancer", "non_cancer")),
                        factor(actual_labels, levels = c("cancer", "non_cancer")))
  
  # Extract performance metrics from confusion matrix
  accuracy <- cm$overall["Accuracy"]
  sensitivity <- cm$byClass["Sensitivity"]
  specificity <- cm$byClass["Specificity"]
  precision <- cm$byClass["Precision"]
  f1 <- cm$byClass["F1"]
  recall <- sensitivity
  
  # Store results for this rule
  rule_results <- rbind(rule_results, data.frame(
    Rule = rule_name,
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1 = f1,
    Recall = recall
  ))
}

# Print the summarized results
print(rule_results)

##Function to plot Confusion Matrices
plot_confusion_matrix <- function(confusion_matrix, title) {
  cm_table <- as.data.frame(confusion_matrix$table)
  ggplot(data = cm_table, aes(x = Prediction, y = Reference)) +
    geom_tile(aes(fill = Freq), color = "white", color = "black") +
    geom_text(aes(label = Freq), fontface = "bold", color = "black", size = 12, vjust = 0.5) +  # Increase text size and adjust positioning
    scale_fill_gradient(low = "lavender", high = "lightblue") +
    labs(title = title, x = "Predicted Label", y = "True Label") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 6, hjust = 0.5),  # Bold and center the title
      axis.title.x = element_text(face = "bold", size = 24, color = "black"),  # Bold x-axis title
      axis.title.y = element_text(face = "bold", size = 24, color = "black"),  # Bold y-axis title
      axis.text.x = element_text(face = "bold", size = 18, color = "black"),   # Bold x-axis text
      axis.text.y = element_text(face = "bold", size = 22, color = "black"),   # Bold y-axis text
      strip.text = element_text(face = "bold", size = 28, color = "black"),    # Bold facet titles
      legend.text = element_text(face = "bold", size = 18),  # Bold and increase legend text size
      legend.title = element_text(face = "bold", size = 16)  # Bold and centered plot title, larger size
    )
}


# 2(E-H) Confusion Matrices for Individual Rules ----------------------


# Create an empty list to store confusion matrices
confusion_matrices <- list()

# Iterate through each rule to compute confusion matrix
for (i in 1:nrow(cancer_rules_train)) {
  # Get genes for the current rule
  gene1 <- cancer_rules_train$Gene1[i]
  gene2 <- cancer_rules_train$Gene2[i]
  rule_name <- cancer_rules_train$Rule[i]
  
  # Extract expression data for the two genes
  expr_gene1 <- expr_subset_train[gene1, ]
  expr_gene2 <- expr_subset_train[gene2, ]
  
  # Apply the rule (e.g., if Gene1 > Gene2, predict "cancer", else predict "no")
  predictions <- as.factor(ifelse(expr_gene1 > expr_gene2, "cancer", "non_cancer"))
  actuals <- as.factor(actual_labels)  # Assuming actual_labels is already a factor
  
  # Calculate the confusion matrix for the current rule
  cm <- confusionMatrix(predictions, actuals)
  
  # Store the confusion matrix in the list
  confusion_matrices[[rule_name]] <- cm
}

# Create a list of confusion matrix plots
cm_plots <- lapply(names(confusion_matrices), function(rule_name) {
  plot_confusion_matrix(confusion_matrices[[rule_name]], rule_name)
})

# Arrange the confusion matrices in a grid
grid <- grid.arrange(grobs = cm_plots, ncol = 2)  # Adjust `ncol` based on the number of rules

#ggsave("Figures/2E_train_bold_confusionmatrices.png", grid, dpi = 600, bg = "white", width = 15, height = 10)

## Expression of each miRNAs between two classes
expr_long_train <- expr_subset_train %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  gather(key = "Sample", value = "Expression", -Gene)

custom_gene_order <- c("4258", "4730", "663a", "1228_5p", "6787_5p", "6800_5p", "1469", "1268a")

# Function for creating BeeswarmPlots
create_beeswarm_plots <- function(dataset, dataset_name) {
  # Initialize list to store beeswarm plots
  beeswarm_plots <- list()
  
  # Loop through each gene in the custom gene order
  for (i in seq_along(custom_gene_order)) {
    gene <- custom_gene_order[i]
    
    # Filter data for the current gene
    gene_data <- dataset %>%
      filter(Gene == gene)
    
    # Ensure gene_data has actual_label and Expression before plotting
    if (nrow(gene_data) > 0) {
      # Create the beeswarm plot for the current gene
      beeswarm_plot <- ggplot(gene_data, aes(x = actual_labels, y = Expression, color = actual_labels)) +
        geom_beeswarm(cex = 0.5, size = 2, alpha = 0.6) +
        facet_wrap(~ Gene, scales = "free_y") +
        scale_color_manual(values = c("cancer" = "red", "non_cancer" = "blue")) +
        labs(x = "Actual Class",
             y = "Expression Level"
        ) +
        theme_minimal(base_size = 8) +
        theme(
          plot.title = element_text(color = "black", face = "bold", hjust = 0.5),
          axis.title.x = element_text(color = "black", face = "bold", size = 8),
          axis.title.y = element_text(color = "black", face = "bold",size = 8),
          axis.text = element_text(color = "black", face = "bold", size = 4),
          strip.text = element_text(color = "black", face = "bold", size = 8),  # Increase the font size of facet labels
          legend.text = element_text(color = "black", face = "bold"),
          legend.title = element_text(color = "black", face = "bold"),
          legend.position = "none"
        )
      
      # Add the plot to the list
      beeswarm_plots[[i]] <- beeswarm_plot
    } else {
      message(paste("No data available for gene:", gene, "in dataset:", dataset_name))
    }
  }
  
  # Check if we have plots to arrange and save
  if (length(beeswarm_plots) > 0) {
    # Create a grid layout (2 rows, 4 columns)
    grid_plot <- grid.arrange(grobs = beeswarm_plots, ncol = 4, nrow = 2)
    
    # Ensure the "train/" directory exists
    dir.create("train", showWarnings = FALSE)
    
    # Save the arranged plots as a single image (grid) in the "train/" directory
    ggsave(filename = paste0("Figures/3ABeeswarm_Plots_", dataset_name, ".png"), 
           plot = grid_plot, width = 24, height = 12, units = "cm")
  }
}



# 3A_Beeswarm plot for miRNAs in training data ------------------------

create_beeswarm_plots(expr_long_train, "Traindata")


#Performance_Summary for Training Data

# Initialize the performance summary data frame
performance_summary <- data.frame(
  Dataset = character(),
  Accuracy = numeric(),
  Sensitivity = numeric(),
  Specificity = numeric(),
  Precision = numeric(),
  F1 = numeric(),
  Recall = numeric(),
  Rule = character(),
  stringsAsFactors = FALSE
)

datasets <- "train"

# Loop through each rule
for (i in 1:nrow(cancer_rules_train)) {
  rule_name <- cancer_rules_train$Rule[i]
  
  # Loop through each dataset
  for (dataset in datasets) {
    
    # Dynamically assign expression data and metadata for the current dataset
    expr_subset <- get(paste0("expr_subset_", dataset))  # expression data
    meta_data <- get(paste0("meta_", dataset))     # metadata
    
    # Assuming 'class' is the column name for actual labels
    actual_labels <- meta_data$class 
    
    # Extract genes for the current rule
    gene1 <- cancer_rules_train$Gene1[i]
    gene2 <- cancer_rules_train$Gene2[i]
    
    # Extract expression data for the two genes
    expr_gene1 <- expr_subset_train[gene1, ]
    expr_gene2 <- expr_subset_train[gene2, ]
    
    # Apply the rule (e.g., if expr_gene1 > expr_gene2, predict "cancer", else predict "no")
    predictions <- ifelse(expr_gene1 > expr_gene2, "cancer", "non_cancer")
    
    # Create confusion matrix
    cm <- confusionMatrix(factor(predictions, levels = c("cancer", "non_cancer")),
                          factor(actual_labels, levels = c("cancer", "non_cancer")))
    
    # Extract performance metrics from confusion matrix
    accuracy <- cm$overall["Accuracy"]
    sensitivity <- cm$byClass["Sensitivity"]
    specificity <- cm$byClass["Specificity"]
    precision <- cm$byClass["Precision"]
    f1 <- cm$byClass["F1"]
    recall <- sensitivity  # Recall is the same as Sensitivity
    
    # Store results for this rule and dataset
    performance_summary <- rbind(performance_summary, data.frame(
      Dataset = dataset,
      Accuracy = accuracy,
      Sensitivity = sensitivity,
      Specificity = specificity,
      Precision = precision,
      F1 = f1,
      Recall = recall,
      Rule = rule_name
    ))
  }
}

# Print the performance summary
print(performance_summary)


# 3B_Performance_Metrices of Individual Rules -------------------------

performance_long <- reshape(performance_summary, 
                            varying = c("Accuracy", "Sensitivity", "Specificity", "Precision", "F1", "Recall"),
                            v.names = "Metric",
                            timevar = "MetricType",
                            times = c("Accuracy", "Sensitivity", "Specificity", "Precision", "F1", "Recall"),
                            direction = "long")


# Convert Rule to a factor with the desired order
performance_long$Rule <- factor(performance_long$Rule, 
                                levels = c("4730,4258", "663a,1228_5p", "6800_5p,6787_5p", "1469,1268a"))

# Plot the heatmap with bold axis labels and black font
ggplot(performance_long, aes(x = Rule, y = MetricType, fill = Metric)) +
  geom_tile() +
  geom_text(aes(label = round(Metric, 3)), color = "black", fontface = "bold") +  # Make text bold
  scale_fill_gradient(low = "lavender", high = "lightblue") +
  labs(title = "Performance Metrics Heatmap for Each Rule", x = "Rule", y = "Metric") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 9, color = "black"),  # Make x-axis labels bold and black
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),  # Make y-axis labels bold and black
    axis.title.x = element_text(face = "bold", size = 14, color = "black"),  # Make x-axis title bold and black
    axis.title.y = element_text(face = "bold", size = 14, color = "black"),  # Make y-axis title bold and black
    strip.text = element_text(face = "bold", size = 12)                     # Make facet strip text bold
  )

ggsave("Figures/3B_performance_metrics_train_individual_ordered.png",height = 6, width = 8)


# 3D_Confusion_matrix_train_for classifier ----------------------------
confusion_train <- caret::confusionMatrix(
  data = factor(results_train$max_score, levels = unique(train_object$data$Labels)),
  reference = factor(meta_train[,"class"], levels = unique(train_object$data$Labels)),
  mode = "everything"
)
print(confusion_train)
plot_confusion_matrix(confusion_train, "Confusion Matrix for Training data")

ggsave("Figures/3D_train_confusion.png")


#  3E Visualize_Classifier_results_in comparison to actual lab --------


platform_colors <- c(
  "benign bone and soft tissue" = "#A3A500",
  "benign breast" = "#A3A500",
  "benign ovary" = "#A3A500",
  "benign prostate" = "#A3A500",
  "biliary tract cancer" = "#E76BF3",
  "bladder cancer" = "#FF61C3",
  "bone and soft tissue sarcomas" = "#8491B4",
  "breast cancer" = "#4E84C4",
  "colorectal cancer" = "#D55E00",
  "esophageal squamous cell cancer" = "#0072B2",
  "extraparenchymal brain tumor and benign brain" = "#CC79A7",
  "gastric cancer" = "#F0E442",
  "hepatocellular cancer" = "#56B4E9",
  "intraparenchymal brain tumors" = "#009E73",
  "lung cancer" = "#E69F00",
  "no cancer" = "#999999",
  "ovarian cancer" = "#D73027",
  "pancreatic cancer" = "#91CF60",
  "prostate cancer" = "#FFBF00"
)

ref_colors <- c("cancer" = "red", "non_cancer" = "blue")

# Plot the heatmap with custom reference label colors
plot_binary_TSP(Data = train_object,        
                classifier = classifier_train,
                binary_col = c("salmon", "lightgreen", "gray"),
                prediction = results_train,
                platform = metadata.subset_GSE211692$state,
                platform_col = platform_colors,
                show_rule_name = TRUE,
                legend = FALSE,
                anno_height = 0.04,
                score_height = 0.075,
                title = "cancer_classifier_train",
                ref_col = ref_colors,
                pred_col = ref_colors,
                margin = c(0,6,0,6)) 

# Create the platform color legend
platform_legend <- Legend(
  title = "Platform/Study",
  at = names(platform_colors),
  legend_gp = gpar(fill = platform_colors),
  nrow = 4,
  title_gp = gpar(fontsize = 14, fontface = "bold", col = "black"),
  labels_gp = gpar(fontsize = 12, fontface = "bold", col = "black")
)
draw(platform_legend)
