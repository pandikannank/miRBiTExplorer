
# 01. Load all required packages ------------------------------------------
setwd("../data/")

##Load all the packages in the environment. 
#If not present, Download from R Base or Bioconductor, whichever is applicable
library(limma)
library(GEOquery) #for downloading GEO datasets
library(umap) #distribution
library(ggplot2)#For visualization
library(ggrepel)#For visualization
library(tidyr)#Operations on dataframes
library(dplyr)#Operations on dataframes
library(multiclassPairs)#SSC
library(data.table)#Operations on dataframes
library(Rtsne)#For visualization
library(caret)#For ML tools_performance Metrices
library(gridExtra)#For Visualization
library(ggbeeswarm)#For Visualization
library(miRNAmeConverter)#miRBase Name conversion
library(RColorBrewer)#For Visualization
library(scales)#For percentage conversion of metrices
library(pROC)#For ROC analysis
library(tibble)#Operations on dataframes
library(Biobase)#Data storage
library(ComplexHeatmap)#For Visualization



# #02_Train ---------------------------------------------------------------

#This code will reproduce the pre-processing of Training Data and Training the classifier
#Fig-2A,2C,2D,2E,2F,2G,2H,3A,3B,3D,3E
########################Preprocessing###########################################

setwd("../data/")

#Training_data_GSE211692
data <- fread("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE211692&format=file&file=GSE211692%5Fprocessed%5Fdata%5Fmatrix%2Etxt%2Egz")
training_data <- as.data.frame(data)
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
metadata.subset_GSE211692 = dplyr::select(metadata_GSE211692,c(2,38))
frequency_table_GSE211692 <- table(metadata.subset_GSE211692$`disease state:ch1`)


# 2A_Disease state Visualization ------------------------------------------
meta_train <- read.csv("train/meta_train_full.csv", row.names = 1, stringsAsFactors = TRUE)
metadisttrain <- metadata.subset_GSE211692
metadisttrain$class <- meta_train$class
metadisttrain$title <- rownames(meta_train)
colnames(metadisttrain)[2] <- "state"

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
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 9, face = "bold", color = "black"),
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

ggsave("../results/Figures/2ASample_distribution_plot.png",height = 10, width = 16)

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

umap_result <- umap(t(highframe_train))# Plot UMAP with specific colors for cancer and non_cancer
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

ggsave("../results/Figures/2Ct-sne_train.png", dpi = 600)


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

ggsave("../results/Figures/2D_mirna_pairs.png", height = 12, width = 12)


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

ggsave("../results/Figures/2E_train_bold_confusionmatrices.png", grid, dpi = 600, bg = "white", width = 15, height = 10)

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
    ggsave(filename = paste0("../results/Figures/3ABeeswarm_Plots_", dataset_name, ".png"), 
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

ggsave("../results/Figures/3B_performance_metrics_train_individual_ordered.png",height = 6, width = 8)


# 3D_Confusion_matrix_train_for classifier ----------------------------
confusion_train <- caret::confusionMatrix(
  data = factor(results_train$max_score, levels = unique(train_object$data$Labels)),
  reference = factor(meta_train[,"class"], levels = unique(train_object$data$Labels)),
  mode = "everything"
)
print(confusion_train)
plot_confusion_matrix(confusion_train, "Confusion Matrix for Training data")

ggsave("../results/Figures/3D_train_confusion.png", dpi = 600)


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

png('../results/Figures/3Efinal4traincomp.png', res = 600, units = "in", width = 16, height = 9, bg = "white")

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

dev.off()

# Create the state color legend
#platform_legend <- Legend(
  #title = "Platform/Study",
  #at = names(platform_colors),
  #legend_gp = gpar(fill = platform_colors),
  #nrow = 4,
  #title_gp = gpar(fontsize = 14, fontface = "bold", col = "black"),
  #labels_gp = gpar(fontsize = 12, fontface = "bold", col = "black")
#)
#draw(platform_legend)




# 03_Test -----------------------------------------------------------------

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
#ggsave("Figures/t-sne_GSE106817.png")

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

#ggsave("Figures/t-sne_GSE112264.png")

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

#ggsave("Figures/t-sne_GSE113486.png")


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

#ggsave("Figures/t-sne_GSE113740.png", dpi = 600, bg = "white")


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

#ggsave("Figures/confusion_GSE113740.png", dpi = 600, bg = "white")

#write.csv(as.data.frame(confusion_GSE113740$table), "final4_confusion_matrix_GSE113740.csv")
#write.csv(meta_GSE113740, "meta_all_GSE113740.csv")

expr_subset_GSE113740 <- GSE113740_frame[rownames(GSE113740_frame) %in% miRNAs_to_plot, ]

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

#ggsave("Figures_t-sne_GSE122497.png")

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
#ggsave("Figures/t-sne_GSE137140.png")


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

#ggsave("Figures/t-sne_GSE139031.png")


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

#ggsave("Figures/t-sne_GSE164174.png")

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
#ggsave("Figures/tsne_GSE73002_nona.png")

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

ggsave("../results/Figures/4A_distribution_test.png", dpi = 600, bg = "white", width = 10, height = 12)


# beeswarm_test -----------------------------------------------------------


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
    ggsave(paste0("../results/Figures/Beeswarm_Plots_", dataset_name, ".png"), grid_plot,
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


ggsave("../results/Figures/4B_Confusion_Matrix_test.png", dpi = 600, bg = "white", height = 12, width = 16)

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

#ggsave("Figures/confusion_merged_Test.png")

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

png('../results/Figures/4C_test_merged.png', res = 600, units = "in", width = 22, height = 10, bg = "white")

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

#platform_legend_test <- Legend(
#title = "Platform/Study",
#at = names(colors),
#legend_gp = gpar(fill = colors),   # Colors for the legend
#ncol = 2,                                  # Number of rows in the legend
#title_gp = gpar(fontsize = 14, fontface = "bold", col = "black"),  # Bold, black title
#labels_gp = gpar(fontsize = 12, fontface = "bold", col = "black")  # Bold, black labels
#)

# Save the platform legend as a separate plot
#png('Figures/test_merged_platform_legend.png', res = 600, units = "in", width = 14, height = 6, bg = "white")
#draw(platform_legend_test)

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


# Reshape the data to long format (if not already done)
cm_metrics_long <- pivot_longer(cm_metrics, cols = -Dataset, names_to = "Metric", values_to = "Value")

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

#ggsave("../results/Figures/4D_heatmap_metrices_bold_ordered.png", dpi = 600, bg = "white", height = 6, width = 8)

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


ggsave("../results/Figures/4D_wb_heatmap_metrices_bold_short.png", dpi = 600, height = 8, width = 10)


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

ggsave("../results/Figures/4E_roc_plots_bold_benign.png", dpi = 600, bg = "white", height = 8, width = 8)




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
  geom_text(aes(label = Freq), color = "black", fontface = "bold", size = 6) +  # Bold the text in tiles
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

ggsave("../results/Figures/5B_Confusion_Matrix_withValidation.png", dpi = 600, height = 8, width = 10)

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



# 05_Benignascancer -------------------------------------------------------

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

#plot_confusion_matrix(confusion_GSE112264_neo, "Confusion Matrix for GSE112264_neo")

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

#ggsave("neoplasm/confusion_GSE113740_neo.png")

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
#ggsave("neoplasm/confusion_GSE137140_neo.png")

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
#ggsave("neoplasm/confusion_GSE139031_neo.png")

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
  geom_text(aes(label = Freq), color = "black", fontface = "bold", size = 6) +  # Bold the text in tiles
  facet_wrap(~ Dataset, ncol = 3) +
  scale_fill_gradient(low = "lavender", high = "lightblue") +
  theme_minimal() +
  labs(title = "Confusion Matrices - Test Datasets",
       x = "Predicted Label", y = "True Label") +
  theme(
    plot.title = element_text(face = "bold", size = 6, hjust = 0.5),  # Bold and center the title
    axis.title.x = element_text(face = "bold", size = 10, color = "black"),  # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 16, color = "black"),  # Bold y-axis title
    axis.text.x = element_text(face = "bold", size = 18, color = "black"),   # Bold x-axis text
    axis.text.y = element_text(face = "bold", size = 18, color = "black"),   # Bold y-axis text
    strip.text = element_text(face = "bold", size = 12, color = "black"),    # Bold facet titles
    legend.text = element_text(face = "bold", size = 18),  # Bold and increase legend text size
    legend.title = element_text(face = "bold", size = 16)  # Bold and increase legend title size
  )

ggsave("../results/Figures/6A_Confusion_benignasnon-cancer_test.png", dpi = 600, bg = "white", height = 12, width = 16)

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
  geom_text(aes(label = Freq), color = "black", fontface = "bold", size = 6) +  # Bold the text in tiles
  facet_wrap(~ Dataset, ncol = 3) +
  scale_fill_gradient(low = "lavender", high = "lightblue") +
  theme_minimal() +
  labs(title = "Confusion Matrices - Test Datasets with Benign as cancer",
       x = "Predicted Label", y = "True Label") +
  theme(
    plot.title = element_text(face = "bold", size = 4, hjust = 0.5),  # Bold and center the title
    axis.title.x = element_text(face = "bold", size = 20, color = "black"),  # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 20, color = "black"),  # Bold y-axis title
    axis.text.x = element_text(face = "bold", size = 5, color = "black"),   # Bold x-axis text
    axis.text.y = element_text(face = "bold", size = 8, color = "black"),   # Bold y-axis text
    strip.text = element_text(face = "bold", size = 12, color = "black"),    # Bold facet titles
    legend.text = element_text(face = "bold", size = 18),  # Bold and increase legend text size
    legend.title = element_text(face = "bold", size = 26)  # Bold and increase legend title size
  )

ggsave("../results/Figures/6B_Confusion_Matrix_withbenignascancer_testlabelsdata.png", dpi = 600, bg = "white", height = 12, width = 14)

# 6C_Train_benignas_non_cancer --------------------------------------------

plot_confusion_matrix(confusion_train, "Confusion Matrix for Training data")
ggsave("../results/Figures/6C_Confusion_benignasnon-cancer_train.png", dpi = 600, height = 8, width = 8)


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

#ggsave("Figures/t-sne_train_with_benign.png")


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
  geom_text(aes(label = Freq), color = "black", fontface = "bold", size = 12) + # Display frequency in tiles
  scale_fill_gradient(low = "lavender", high = "lightblue") +
  theme_minimal() +
  labs(title = "Confusion Matrix for Neo Datasets",
       x = "Predicted Label", y = "True Label") +
  theme(
    plot.title = element_text(face = "bold", size = 8, hjust = 0.5),  # Bold and center the title
    axis.title.x = element_text(face = "bold", size = 16, color = "black"),  # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 16, color = "black"),  # Bold y-axis title
    axis.text.x = element_text(face = "bold", size = 14, color = "black"),   # Bold x-axis text
    axis.text.y = element_text(face = "bold", size = 14, color = "black"),   # Bold y-axis text
    strip.text = element_text(face = "bold", size = 20, color = "black"),    # Bold facet titles
    legend.text = element_text(face = "bold", size = 16),  # Bold and increase legend text size
    legend.title = element_text(face = "bold", size = 16)  # Bold and increase legend title size
  )
ggsave("../results/Figures/6D_merged_confusion_Train_benign as cancer.png", dpi = 600, height = 12, width = 12)


# 6E_Test_Benignasnon_cancer_merged_conf ----------------------------------

plot_confusion_matrix(confusion_test_merged, "Confusion Matrix for Test_Datasets")
ggsave("../results/Figures/6E_confusion_merged_Test.png", dpi = 600, height = 12, width = 12)


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

ggsave("../results/Figures/6F_merged_confusion_Test_benign as cancer.png", dpi = 600, height = 12, width = 12)

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

ggsave("../results/Figures/6G_heatmap_metrices_benignneo.png",, dpi = 600, height = 10, width = 20)


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

ggsave("../results/Figures/6H_roc_plots_bold_benignascancer.png", dpi = 600, bg = "white", height = 8, width = 16)
