# pcitRif

## Installation

<p>To properly run <b>pcitRif</b> package is necessary to install some dependencies:</p>
<p>for Linux users is necessary to install the dependencies:</p>

* libcurl4-openssl-dev
* libxml2-dev 
* libssl-dev

<p>To install R packages dependencies, run:</p>

```R
#CRAN dependencies
packagesCRAN <- c("crayon", "pbapply", "reshape2", "kableExtra", "knitr", "rmarkdown", "ggplot2", "gridExtra", "BiocManager")
install.packages(packagesCRAN[!packagesCRAN %in% installed.packages()[,1]])

#Bioconductor dependencies
packagesBioc <- c("airway", "SummarizedExperiment", "DESeq2")
BiocManager::install(packagesBioc[!packagesBioc %in% installed.packages()[,1]])
```
