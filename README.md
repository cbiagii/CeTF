# pcitRif <img src="vignettes/logo.png" align="right" width="300" />
An implementation of PCIT and RIF analysis in R



## Overview
This package provides the necessary instructions for performing the Partial Correlation coefficient with Information Theory (PCIT) from [Reverter and Chan 2008](https://doi.org/10.1093/bioinformatics/btn482) and Regulatory Impact Factors (RIF) from [Reverter et al. 2010](https://doi.org/10.1093/bioinformatics/btq051) algorithm. The PCIT algorithm identifies meaningful correlations to define edges in a weighted network. The algorithm can be applied to any correlation-based network including but not limited to gene co-expression networks. While the RIF algorithm identify critical transcript factors (TF) from gene expression data. These two algorithms when combined provide a very relevant layer of information for gene expression studies (Microarray, RNA-seq and single-cell RNA-seq data).

## Installation

There are three ways to install this package, conventional, docker and using BioConductor installation.

### Conventional

To properly run <b>pcitRif</b> package is necessary to install some dependencies:
for Linux users is necessary to install the dependencies:

* libcurl4-openssl-dev
* libxml2-dev 
* libssl-dev
* gfortran
* build-essential
* libz-dev 
* zlib1g-dev


To install R packages dependencies, run:

```R
#CRAN dependencies
packagesCRAN <- c("crayon", "pbapply", "reshape2", "kableExtra", "knitr", "rmarkdown", "ggplot2", "gridExtra", "BiocManager")
install.packages(packagesCRAN[!packagesCRAN %in% installed.packages()[,1]])

#Bioconductor dependencies
packagesBioc <- c("airway", "SummarizedExperiment", "DESeq2")
BiocManager::install(packagesBioc[!packagesBioc %in% installed.packages()[,1]])

```

Finally, to install pcitRif package:

```R
devtools::install_github("cbiagii/pcitRif")
```

### Docker

To install docker follow the instructions in the links below depending on your operating system

* CentOS: https://docs.docker.com/install/linux/docker-ce/centos/
* Debian: https://docs.docker.com/install/linux/docker-ce/debian/
* Fedora: https://docs.docker.com/install/linux/docker-ce/fedora/
* Ubuntu: https://docs.docker.com/install/linux/docker-ce/ubuntu/
* MacOS: https://docs.docker.com/docker-for-mac/install/
* Windows: https://docs.docker.com/docker-for-windows/install/

#### Docker pull

Once docker is installed, the next step is pull the target-prediction image from dockerhub using the following command:

```docker
docker pull biagii/pcitRif
```

#### Running image

There are several differents parameters to run the downloaded image. The most commom way is executing the following command:

```docker
docker run --rm --name [ANY_NAME] -v /server/path/:/docker/path/ biagii/pcitRif
```
## Help

<p>Any questions contact the developer by email: <a href="#">biagi@usp.br</a></p>
