# miRBiTExplorer
This repository contains codes and files to reproduce the analyses performed in the Manuscript titled, "Global Pan-cancer serum miRNA classifier across 13 cancer types: Analysis of 46,349 clinical samples".
Please run all the codes in the sequence.
Before starting, download all the files from https://zenodo.org/records/14906312 (Publicly available after Publication) and unzip it.

01_Packages.R contains all the necessary R libraries that needs to be loaded to perform the analysis. These packages can be downloaded from https://stat.ethz.ch/R-manual/R-devel/library/base/html/00Index.html or Bioconductor https://www.bioconductor.org/
02_Train.R contains the code to reproduce the pre-processing of Training Data and single-sample classifier training. Fig-2A,2C,2D,2E,2F,2G,2H,3A,3B,3D,3E can be generate with the code.
03_Test.R contains the code to reproduce the pre-processing of Test Data and the classifier performance. Fig-4A,4B,4C,4D,4E can be generated with this code.
04_Validation.R contains the code to reproduce the reproduce the pre-processing of Validation Data and the classifier performance. Fig-5A,5B,5D,5F,5H can be generated with this code.
05_Benign_as_cancer.R contains the code to reproduce the analysis when benign is considered cancer/benign class and comparisons between cancer and cancer/benign groups. Fig-6A,6B,6C,6D,6E,6F,6G,6H can be generated with this code.
Fig.1, Fig.2B, Fig. 3C, Fig. 5C, Fig.5E, Fig.5G are schematics or made with simple tools and hence their codes are not applicable to be present here.
