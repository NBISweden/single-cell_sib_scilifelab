---
title: "Batch correction for scRNA-seq data"
author: "Panagiotis Papasaikas"
date: "2019-10-16"
output:
  rmarkdown::html_document:
    highlight: pygments
    toc: true
    toc_depth: 3
    keep_md: yes
editor_options: 
  chunk_output_type: console
bibliography: batch_correction.bib
---




# Introduction

In this lab we will focus on  data integration / batch correction apporaches
specifically appropriate for single cell RNAseq datasets. 
We will go through the steps of 
1. batch effect diagnosis, 
2. actual correction
3. evaluation of the effects/quality of correction.

The dataset that we will use is a composite dataset of three independent
10x runs originating from different labs. It consists of 9288 mammary 
epithelial cells,sequenced using 10x Genomics technology, which has already
been pre-filtered to include only cells that are assigned unambiguously
to one of three major  cell types:
luminal progenitors, luminal mature and basal.


covers some of the most commonly used methods for finding
differentially expressed genes ("marker genes") between clusters in single-cell
RNA-seq. We will use an example data set consisting of 2,700 PBMCs, sequenced
using 10x Genomics technology.

Several of the methods outlined here are based on the batchelor package developed
by the lab of Aaron lun and thus the corresponding R vignette and manual are excellent
sources of additional information ["Single-cell correction with batchelor"](https://bioconductor.org/packages/release/bioc/vignettes/batchelor/inst/doc/correction.html).
The student can also refer to the Integrating datasets chapther of ["Orchestrating single-cell analysis with Bioconductor"](https://osca.bioconductor.org/integrating-datasets.html).
which is also largely based on the batchelor package.

# Load packages

We first load the required R packages. 


```r
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(scater)
  library(ggplot2)
  library(scran)
  library(batchelor)
  library(BiocSingular)
  library(BiocNeighbors)
  library(e1071)
  library(coop)
  library(Rtsne)
  source("Helper_Functions.R")
})
```

# Load the pre-processed dataset:

Next, we load the preprocessed dataset and have a first look at the 
composition of the dataset:


```r
## Download the data and set row names to gene symbols whenever possible
bct <- readRDS(gzcon(url("https://github.com/NBISweden/single-cell_sib_scilifelab/blob/master/datasets/SCE_MammaryEpithelial_x3.rds?raw=true")))
rownames(bct) <- scater::uniquifyFeatureNames(
  ID = rownames(bct), 
  names = as.character(rowData(bct)$gene.symbols)
)

## Dataset compostion per cell type and study:  
table(colData(bct)$study , colData(bct)$cell.class)
```

```
##      
##       luminal_progenitor basal luminal_mature
##   spk               1097   475           1043
##   vis                694  1257           1096
##   wal                729   780           2117
```


# Batch effect diagnosis:

We will now look at some high level characteristics of the three datasets
that clearly indicate the presence of batch effects:




```r
## Check library size amd  number of detected gene distribution in the three datasets:
   ggplot(data.frame(colData(bct)), 
                 aes(x = study, y = library_size, color=study)) + 
              geom_violin() + theme_bw() + 
              ylab("Total UMI counts per cell") + 
              ggtitle("Library size distribution per study") 
```

![](batch_correction_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
   ggplot(data.frame(colData(bct)), 
              aes(x = study, y = detected_genes, color=study)) + 
              geom_violin() + theme_bw() + 
              ylab("NODG") + 
              ggtitle("Number of detected genes per study") 
```

![](batch_correction_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

We can already see the the 3 studies differ in these high level characterstics.
Let's now take a look at the impact of the batch in the TSNE projections generated
without taking any steps for batch correction (apart from a simple library normalization ):



```r
# We first normalize all cells for library size.
assays(bct )[["lognorm"]] <- lin_norm_log(as.matrix(  logcounts(bct) ), factor=1e+04)  
reducedDim(bct, "PCA" )  <- rsvd::rpca(t( bct@assays[["lognorm"]]  ),k=32,retx=TRUE,center=TRUE,scale=FALSE)$x
reducedDim(bct, "TSNE" ) <- Rtsne(bct@reducedDims$PCA, perplexity = 30, initial_dims=32, pca=FALSE, theta=0.3)$Y #~15-60 seconds run time

cowplot::plot_grid(scater::plotTSNE(bct, colour_by = "study" ),
                   scater::plotTSNE(bct, colour_by = "cell.class"))
```

![](batch_correction_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

# Library normalization and PCA methods customized for the presence of batches

We saw above that our data form distinct clusters according for both study origin as well as biological cell type.
We now take a first step towards batch correction using approaches for 
A. library normalization and 
B. PCA 
that explicitly account for the presence of batches.


The new library normalization strategy is essentially a 3-step process:
1. Similar cells (cell with similar expression profiles) are pooled together using the scran function quickCluster (or any othe sensible clustering method).
2. Size factors are computed using the deconvolution method (Lun et al., 2016) and using the scran function computeSumFactors
3. Systematic differences in coverage across batches are removed by rescaling the size factors using median-based
normalization on the ratio of the average counts between batches.This is done using the multiBatchNorm() function  from the batchelor package. 
The function multiBatchNorm from the scran package provides identical functionality.

The new PCA strategy is performed using the multibatchPCA function from the batchelor package. This is essentially
performing PCA on the merged matrxi of the three batches but now each batch is "forced" to contribute equally to the
identification of the loading vectors. In addition each batch's contribution to the gene-gene covariance matrix is 
normalized by the corresponding number of cells.


```r
bct.mBN <- bct #Create a copy of the single cell experiment

# Clustering:
clusters <- scran::quickCluster(bct.mBN , method = "igraph", min.mean = 0.1,
                                  use.ranks = FALSE, BSPARAM = IrlbaParam(),
                                  irlba.args = list(maxit = 1000),
                                  #BPPARAM = MulticoreParam(workers = 16)
                                 )
print(table(clusters))
```

```
## clusters
##    1    2    3    4    5    6    7    8    9   10   11   12   13 
## 1112 1032 1253  608  849  692 1280  794  331  707  200  204  226
```

```r
# Size factor calculation using the deconvolution method:
bct.mBN <- scran::computeSumFactors(bct.mBN , min.mean = 0.1, cluster = clusters,
                                #BPPARAM = MulticoreParam(workers = 16)
                                ) # 30"-60" running time


# Library normalization that corrects for batch-specific difference in size factors:
mBN <- batchelor::multiBatchNorm(V = bct.mBN[,colData(bct.mBN )$study=="vis"],
                                 W = bct.mBN[,colData(bct.mBN )$study=="wal"],
                                 S = bct.mBN[,colData(bct.mBN )$study=="spk"])
bct.mBN <- cbind( mBN$V, mBN$W, mBN$S) 
# Notice that only the logcount slot that is affected. 
# The counts (and the custom lognorm assays) of the original bct object 
# are unaffected.

# Mutli-batch PCA that makes sure each batch contributes equally to the loading vectors:
mB.PCA <- batchelor::multiBatchPCA( bct.mBN, batch=colData(bct.mBN)$study, d=32, preserve.single = TRUE)
```

```
## Warning in sweep(centered, 2, w, "/", check.margin = FALSE): 'check.margin' is ignored when 'x' is a DelayedArray object or
##   derivative
```

```r
reducedDim(bct.mBN , "PCA" )  <- mB.PCA[[1]]

# Let's now recalculate the projection and see if the improved normalization
# affected the quality of the projections:
reducedDim(bct.mBN , "TSNE" ) <- Rtsne(bct.mBN@reducedDims$PCA, perplexity = 30, initial_dims=32, pca=FALSE, theta=0.3)$Y #~15-60 seconds run time

cowplot::plot_grid(scater::plotTSNE(bct.mBN, colour_by = "study" ),
                   scater::plotTSNE(bct.mBN, colour_by = "cell.class")
                   )
```

![](batch_correction_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

# Batch effect correction using linear regression

The most common approaches for correcting batch effects in the case of bulk RNA-seq datasets are based on linear regression.
Typically a linear model is fitted per-gene where the batch effect is accounted for by a separate term in our model specification.
One we fit thie model, the batch-specfic effects can be corrected by setting this term to zero.
This is the basis of the removeBatchEffect() function from the limma package (Ritchie et al. 2015) as well the comBat() function
from the sva package (Leek et al. 2012).
The rescaleBatches function from the batchelor package takes a similar approach, with adjustment for improving efficiency in the
single-cell regime. Specifically, for each gene, the mean expression in each batch is scaled down until it is equal to the 
lowest mean across all batches. 



```r
bct.linCor <- bct.mBN 
bct.linCor  <- batchelor::rescaleBatches(bct.linCor, batch=colData(bct.linCor)$study)
colData(bct.linCor) <-  colData(bct.mBN )  


reducedDim(bct.linCor, "PCA" )  <- rsvd::rpca(t( bct.linCor@assays[["corrected"]]  ),k=32,retx=TRUE,center=TRUE,scale=FALSE)$x
reducedDim(bct.linCor, "TSNE" ) <- Rtsne(bct.linCor@reducedDims$PCA, perplexity = 30, initial_dims=32, pca=FALSE, theta=0.3)$Y #~15-60 seconds run time

cowplot::plot_grid(scater::plotTSNE(bct.linCor, colour_by = "study" ),
                   scater::plotTSNE(bct.linCor, colour_by = "cell.class")
                   )
```

![](batch_correction_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

An interesting exercise is to repeat the linear correction, but this time using as input the bct object where
library normalization has not taken into account the presence of batches. Do you see any differences?



# Batch effect correction using fast-MNN




```r
d <- 32
FMNN.out <-  batchelor::fastMNN( bct.mBN  , batch=bct.mBN$study , use.dimred="PCA", d=d ) 
```

```
## Warning in (function (jobs, data, centers, info, distance, k, query,
## get.index, : tied distances detected in nearest-neighbor calculation
```

```r
reducedDim (bct.mBN, "PCA.FMNN" ) <- FMNN.out$corrected 

reducedDim(bct.mBN, "TSNE" ) <- Rtsne( bct.mBN@reducedDims$PCA.FMNN, perplexity = 30, initial_dims=64, pca=FALSE,num_threads=32,theta=0.25)$Y

cowplot::plot_grid(scater::plotTSNE(bct.mBN, colour_by = "study" ),
                   scater::plotTSNE(bct.mBN, colour_by = "cell.class")
                   )
```

![](batch_correction_files/figure-html/unnamed-chunk-7-1.png)<!-- -->



# Assessment of the quality of batch correction



# Session info


```r
sessionInfo()
```

```
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS Mojave 10.14.5
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] Rtsne_0.15                  coop_0.6-2                 
##  [3] e1071_1.7-2                 BiocNeighbors_1.2.0        
##  [5] BiocSingular_1.0.0          batchelor_1.0.1            
##  [7] scran_1.12.1                scater_1.12.2              
##  [9] ggplot2_3.2.1               SingleCellExperiment_1.6.0 
## [11] SummarizedExperiment_1.14.1 DelayedArray_0.10.0        
## [13] BiocParallel_1.18.1         matrixStats_0.55.0         
## [15] Biobase_2.44.0              GenomicRanges_1.36.1       
## [17] GenomeInfoDb_1.20.0         IRanges_2.18.3             
## [19] S4Vectors_0.22.1            BiocGenerics_0.30.0        
## [21] BiocStyle_2.12.0           
## 
## loaded via a namespace (and not attached):
##  [1] viridis_0.5.1            dynamicTreeCut_1.63-1   
##  [3] edgeR_3.26.8             viridisLite_0.3.0       
##  [5] DelayedMatrixStats_1.6.1 assertthat_0.2.1        
##  [7] statmod_1.4.32           BiocManager_1.30.7      
##  [9] dqrng_0.2.1              GenomeInfoDbData_1.2.1  
## [11] vipor_0.4.5              yaml_2.2.0              
## [13] pillar_1.4.2             lattice_0.20-38         
## [15] glue_1.3.1               limma_3.40.6            
## [17] digest_0.6.21            XVector_0.24.0          
## [19] colorspace_1.4-1         cowplot_1.0.0           
## [21] htmltools_0.4.0          Matrix_1.2-17           
## [23] pkgconfig_2.0.3          zlibbioc_1.30.0         
## [25] purrr_0.3.2              scales_1.0.0            
## [27] tibble_2.1.3             withr_2.1.2             
## [29] lazyeval_0.2.2           magrittr_1.5            
## [31] crayon_1.3.4             evaluate_0.14           
## [33] class_7.3-15             beeswarm_0.2.3          
## [35] tools_3.6.1              stringr_1.4.0           
## [37] munsell_0.5.0            locfit_1.5-9.1          
## [39] irlba_2.3.3              compiler_3.6.1          
## [41] rsvd_1.0.2               rlang_0.4.0             
## [43] grid_3.6.1               RCurl_1.95-4.12         
## [45] igraph_1.2.4.1           bitops_1.0-6            
## [47] labeling_0.3             rmarkdown_1.16          
## [49] gtable_0.3.0             R6_2.4.0                
## [51] gridExtra_2.3            knitr_1.25              
## [53] dplyr_0.8.3              stringi_1.4.3           
## [55] ggbeeswarm_0.6.0         Rcpp_1.0.2              
## [57] tidyselect_0.2.5         xfun_0.10
```

# References

