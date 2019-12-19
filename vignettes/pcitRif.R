## ----setup, echo=FALSE, results="hide"----------------------------------------
knitr::opts_chunk$set(tidy = FALSE,
                      cache = FALSE,
                      dev = "png",
                      message = FALSE, error = FALSE, warning = TRUE)

## ----install1, eval=FALSE-----------------------------------------------------
#  #CRAN dependencies
#  packagesCRAN <- c("crayon", "geomnet", "GGally", "ggplot2", "ggpubr", "graphics", "kableExtra", "knitr", "network", "pbapply", "reshape2", "rmarkdown", "scales", "stats", "testthat", "utils")
#  install.packages(packagesCRAN[!packagesCRAN %in% installed.packages()[,1]])
#  
#  #Bioconductor dependencies
#  packagesBioc <- c("airway", "clusterProfiler", "DESeq2", "org.Hs.eg.db", "SummarizedExperiment")
#  BiocManager::install(packagesBioc[!packagesBioc %in% installed.packages()[,1]])

## ----install2, eval=FALSE-----------------------------------------------------
#  install.packages("/path/to/package/pcitRif_0.1.0.tar.gz")

## ----PCIT, eval=TRUE, warning=F, message=F------------------------------------
# Loading packages
library(pcitRif)
library(airway)
library(kableExtra)
library(knitr)

# Loading airway data
data("airway")

# Creating a variable with annotation data
anno <- as.data.frame(colData(airway))
anno <- anno[order(anno$dex, decreasing = TRUE), ]
anno <- data.frame(cond = anno$dex, 
                   row.names = rownames(anno))

# Creating a variable with count data
counts <- assay(airway)

# Sorting count data samples by conditions (untrt and trt)
counts <- counts[, rownames(anno)]
colnames(counts) <- paste0(colnames(counts), c(rep("_untrt", 4), rep("_trt", 4)))

# Differential Expression analysis to use only informative genes
DEGenes <- expDiff(exp = counts,
                   anno = anno,
                   conditions = c('untrt', 'trt'),
                   lfc = 4,
                   padj = 0.05, 
                   diffMethod = "Reverter")

# Selecting only DE genes from counts data
counts <- counts[rownames(DEGenes$DE_unique), ]

# Converting count data to TPM
tpm <- countsToTPM(counts)

# Count normalization
PCIT_input <- normExp(tpm)

# PCIT input for untrt
PCIT_input_untrt <- PCIT_input[,grep("_untrt", colnames(PCIT_input))]

# PCIT input for trt
PCIT_input_trt <- PCIT_input[,grep("_trt", colnames(PCIT_input))]

# Performing PCIT analysis for untrt condition
PCIT_out_untrt <- PCIT(PCIT_input_untrt, tolType = "mean")

# Performing PCIT analysis for trt condition
PCIT_out_trt <- PCIT(PCIT_input_trt, tolType = "mean")

# Printing first 10 rows for untrt condition
kable(PCIT_out_untrt$tab[1:10, ]) %>% 
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

# Printing first 10 rows for trt condition
kable(PCIT_out_trt$tab[1:10, ]) %>% 
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

# Adjacency matrix: PCIT_out_untrt$adj_raw or PCIT_out_trt$adj_raw
# Adjacency matrix only with the significant values: PCIT_out_untrt$adj_sig or PCIT_out_trt$adj_sig

## ----hist, eval=TRUE, fig.align="center"--------------------------------------
# Example for trt condition
histPlot(PCIT_out_trt$adj_sig)

## ----densityRaw, eval=TRUE, fig.align="center"--------------------------------
# Example for trt condition
densityPlot(PCIT_out_trt$adj_raw)

## ----densitySig, eval=TRUE, fig.align="center"--------------------------------
# Example for trt condition
densitySig(mat1 = PCIT_out_trt$adj_raw, 
           mat2 = PCIT_out_trt$adj_sig, 
           threshold = 0.5)

## ----RIF, eval=TRUE, warning=F, message=F-------------------------------------
# Loading packages
library(pcitRif)
library(airway)
library(kableExtra)
library(knitr)

# Loading airway data
data("airway")

# Creating a variable with annotation data
anno <- as.data.frame(colData(airway))
anno <- anno[order(anno$dex, decreasing = TRUE), ]
anno <- data.frame(cond = anno$dex, 
                   row.names = rownames(anno))

# Creating a variable with count data
counts <- assay(airway)

# Sorting count data samples by conditions (untrt and trt)
counts <- counts[, rownames(anno)]
colnames(counts) <- paste0(colnames(counts), c(rep("_untrt", 4), rep("_trt", 4)))

# Differential Expression analysis to use only informative genes
DEGenes <- expDiff(exp = counts,
                   anno = anno,
                   conditions = c('untrt', 'trt'),
                   lfc = 1.5,
                   padj = 0.05, 
                   diffMethod = "Reverter")

# Selecting only DE genes from counts data
counts <- counts[rownames(DEGenes$DE_unique), ]

# Converting count data to TPM
tpm <- countsToTPM(counts)

# Count normalization
Clean_Dat <- normExp(tpm)

# Loading the Transcript Factors (TFs) character
data("TFs")

# Verifying which TFs are in the subsetted normalized data
TFs <- rownames(Clean_Dat)[rownames(Clean_Dat) %in% TFs]

# Selecting the Target genes
Target <- setdiff(rownames(Clean_Dat), TFs)

# Ordering rows of normalized count data
RIF_input <- Clean_Dat[c(Target, TFs), ]

# Performing RIF analysis
RIF_out <- RIF(input = RIF_input,
               nta = length(Target),
               ntf = length(TFs),
               nSamples1 = 4,
               nSamples2 = 4)

# Printing first 10 rows
kable(RIF_out[1:10, ]) %>% 
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

## ----all, eval=TRUE-----------------------------------------------------------
# Loading packages
library(airway)
library(pcitRif)

# Loading airway data
data("airway")

# Creating a variable with annotation data
anno <- as.data.frame(colData(airway))

# Creating a variable with count data
counts <- assay(airway)

# Sorting count data samples by conditions (untrt and trt)
counts <- counts[, order(anno$dex, decreasing = TRUE)]

# Loading the Transcript Factors (TFs) character
data("TFs")

# Performing the complete analysis
out <- runAnalysis(mat = counts, 
                   conditions=c("untrt", "trt"),
                   lfc = 3.5,
                   padj = 0.05,
                   TFs = TFs,
                   nSamples1 = 4,
                   nSamples2= 4,
                   tolType = "mean",
                   diffMethod = "Reverter", 
                   data.type = "counts")

## ----SmearPlotDE, eval=TRUE, fig.align="center"-------------------------------
# Using the runAnalysis output (pcitrif class object)
SmearPlotDE(object = out,
            diffMethod = "Reverter",
            lfc = 1.5,
            conditions = c("untrt", "trt"))

## ----SmearPlotTF, eval=TRUE, fig.align="center"-------------------------------
# Using the runAnalysis output (pcitrif class object)
SmearPlotTF(object = pcitrifExample,
            diffMethod = "Reverter",
            lfc = 1.5,
            conditions = c("untrt", "trt"),
            TF = "ENSG00000185917")

## ----singleNetworkPlot, eval=TRUE, fig.align="center"-------------------------
# Using the runAnalysis output (pcitrif class object)
singleNetworkPlot(out)

## ----getGroupGO, eval=TRUE, fig.align="center"--------------------------------
# Loading Homo sapiens annotation package
library(org.Hs.eg.db)

# Accessing the network for condition 1
genes <- unique(c(as.character(getNet1(out)[,1]), 
                  as.character(getNet1(out)[,2])))

# Performing getGroupGO analysis
cond1 <- getGroupGO(genes = genes, 
                    ont = "BP", 
                    keyType = "ENSEMBL", 
                    annoPkg = org.Hs.eg.db)

