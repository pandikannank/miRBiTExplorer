
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