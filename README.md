# StatinEMT


## Introduction
This is a project to regenerate the pharmacogenomics analysis in this paper:
"Statin-induced cancer cell death can be mechanistically uncoupled from prenylation of RAS family proteins"
citation's info to be added

The script for the analysis is written in R.


----

## Dependencies:
These R packages need to be installed in order to run the analysis script:
- Biobase
- PharmacoGx
- mixtools
- mclust


----
## Reproducibility of the Analysis:
- Once the project is downloaded to the user computer, the user needs to navigate to the main directory of the project "StatinEMT-master".
- Inside the main directory, there is an R script file named "Script_EMT_Figures_Informatics.R". Running this script will regenerate the full pharmacogenomics analysis of the paper.

**Important Note:** the user needs to set the working directory inside the script file before running it, i.e. change the following code inside the script:

`StatinEMT_Directory <- NULL`
to
`StatinEMT_Directory <- "./pathToStatinEMT-master_Directory"`

- Running the script will produce all the figures related to the pharmacogenomics analysis
- Included in the project are the bimodality scores for all genes across CCLE cell lines that have RNA-seq data. This can be found in "./data/BimodalScores_All_Genes.RData"

