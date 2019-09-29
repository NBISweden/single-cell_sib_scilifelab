## This file is read by install_packages.R to automatically prepare
## the R environment for the course.
## - please add your requirements below (R packages, additional code to run)
## - please don't modify the lines containing START-PACKAGES, END-PACKAGES,
##   START-EXPRESSIONS, END-EXPRESSIONS

## R packages: List here required packages for exercises for R 3.6.1, Bioc 3.9
## (one package per line in format: pkgname[-minimal-version])
## START-PACKAGES:
AUCell
BiocSingular
CATALYST
cluster
clustree
dendextend
dynamicTreeCut
dyno
edgeR
flowCore
HDCytoData
igraph
iSEE
limma
MAST
mbkmeans
mclust
pheatmap
SC3
scater
scran
Seurat-3.1
SingleCellExperiment
TENxPBMCData
tidyverse
zinbwave
## END-PACKAGES

## Additional code: Add here additional expressions to be run, e.g. to download
## data sets from ExperimentHub (or other sources) in advance
## START-EXPRESSIONS:
TENxPBMCData::TENxPBMCData(dataset = "pbmc3k")
## END-EXPRESSIONS
