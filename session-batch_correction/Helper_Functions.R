# @param m  a (potentially sparse) gene x cells count matrix
# @param f  a number between 0 and 1 the fraction of overdispersed genes to keep
# @return a vector of rownames of genes to keep.
select_variable_genes<-function(m,f) {
  zeroes=which(rowSums(m) <= max(1,min(rowSums(m))) )
  all.nz.genes   <- rownames(m[-zeroes,])
  m <- m[all.nz.genes,]
  df <- data.frame(mean = rowMeans(m + 1/ncol(m)), cv = apply(m,1,sd) / rowMeans(m + 1/ncol(m)), var = apply(m,1,var))
  df$dispersion <- with(df, var/mean)
  df$mean_bin <- with(df, cut(mean, breaks = c(-Inf, unique(quantile(mean, seq(0.1,1,0.02), na.rm = TRUE)), Inf)))
  var_by_bin <- data.frame(mean_bin = factor(levels(df$mean_bin), levels = levels(df$mean_bin)),
                           bin_median = as.numeric(tapply(df$dispersion, df$mean_bin, stats::median)),
                           bin_mad = as.numeric(tapply(df$dispersion, df$mean_bin, stats::mad)))[table(df$mean_bin) > 0,]
  df$bin_disp_median <- var_by_bin$bin_median[match(df$mean_bin, var_by_bin$mean_bin)]
  df$bin_disp_mad <- var_by_bin$bin_mad[match(df$mean_bin, var_by_bin$mean_bin)]
  df$dispersion_norm <- with(df, (dispersion - bin_disp_median)/(bin_disp_mad + 0.01) )
  
  n_genes_keep=ceiling(f*nrow(m) ) #In the end retain only the top 100*f% overdispersed genes
  disp_cut_off <- sort(df$dispersion_norm,decreasing=TRUE)[n_genes_keep]
  genes_keep <- which(df$dispersion_norm >= disp_cut_off)
  return(rownames(m)[genes_keep])
}



# Convert log transformed data back to linear space / library normalize
lin_norm_log <- function(x, library.normalize=TRUE, log=TRUE, factor=1e+04) { 
  x <- 2^x-1
  if (library.normalize){
    x <-   sweep(x, 2, colSums(x), FUN = "/") * factor 
  }
  x[x < 1e-2]=0
  if (log){
    x <- log2(x+1)  
  }
  x
}


# Map value vector to colors
map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}


# Calculates KNN smoothened expression values for one or more genes.
# Uses a reduced dim representation for KNN identification {proj: cells x features }
# and an expression matrix (typically logexpr) for gene expression value lookup {exp.matrix: fetures x cells }
# markers is a vector of indices/rownames corresponding to the rows of the expression matrix
# that contain the genes of interest
# returns a numeric vector of length cells
knn.smooth <-function(proj,exp.matrix,  markers, k=50 ){
library(BiocNeighbors)
NN <- BiocNeighbors::findKmknn( proj, k=k , BPPARAM=MulticoreParam(workers=8) )
if (length(markers)==1){
    marker.vals <- exp.matrix[markers,as.vector(NN$index)]
} else{
    marker.vals  <- colMeans(exp.matrix[markers,])
    marker.vals <- marker.vals[as.vector(NN$index)]
}
marker.vals <- matrix(marker.vals,nrow=nrow(NN$index),byrow = FALSE)
marker.levels <- rowMeans(marker.vals)
}






##### Functions for evaluating bc performance:

#### A. Mixing:
### Check if a simple svm model can discriminate between batches
### Input is a K x d matrix X with K samples and d features and an annotation vector Y
### n: Number of samples to use
### seed: FALSE or a number for setting the random seed
### ... Additional arguments passed to e1071::svm
### Returns cross validation classification accuracy (0-100).
### Lower is better
library(e1071)
mixing.svm <- function (X, Y, n=1000, seed=1, balanced=TRUE, ...) {
    svm.args <- list(...)
    if (is.null(svm.args$kernel)) svm.args$kernel <- "radial" # use kernel="linear" for linear SVM. kernel="radial" for non-linear SVM
    if (is.null(svm.args$cost))   svm.args$cost <- 10
    if (is.null(svm.args$scale))  svm.args$scale <- TRUE
    if (is.null(svm.args$cross))  svm.args$cross <- 5
    if (seed){
        set.seed(seed)
    }
    
    if (balanced) {
        nb <- min( floor( n/length(unique(Y) ) ) , min(table(Y))  )
        s <- sapply(unique(Y), function(x) sample( which(Y==x), nb ),simplify = TRUE   )
    }  else{  
        s <- sample(1:nrow(X), min(n, nrow(X)) )
    }
    
    dat <- data.frame( x=X[s,], y=as.factor(Y[s]) )  
    
    svmfit <- do.call(svm, c( list(y ~ . , data=dat ),svm.args)     )
    accuracy <- svmfit$tot.accuracy
    return(accuracy)
}
#e.g test <- mixing.svm(M3,study_annot,n=700, balanced=TRUE)



#### B1. Structure preservation. Local:
### Input are two K x d matrices M1, M2 with K samples and d features.
### n is the sample size of K
### Average Jaccard distance between a sample of n vertices before and after mixing
### Ranges from 0-1. Lower is better.
library(BiocNeighbors)
library(scran)
local.dist <- function (m1, m2, metric=c("Jaccard"), n=1000, k=20,  seed=1, cosnorm=FALSE ) {
    
    if (!identical(dim(m1), dim(m2))) {
        stop("The two input matrices must have the same dimensions")
    }
    
    if (seed){
        set.seed(seed)
    }
    
    s <- sample(1:nrow(m1), min(n, nrow(m1)) )
    m1  <- m1[s,]
    m2  <- m2[s,] 
    
    # Cosine normalization of feature matrices:
    if (cosnorm) {
        m1 <-   t(scran::cosineNorm( t(m1) ))
        m2 <-   t(scran::cosineNorm( t(m2) ))
    }
    
    #Find kNNs for each point in m1, m2
    N1 <- BiocNeighbors::findKNN(m1, k=k, get.dist=FALSE)$index
    N2 <- BiocNeighbors::findKNN(m2, k=k, get.dist=FALSE)$index
    
    #Calculate average Jaccard similarity:
    Jindex <- mean ( sapply(1:nrow(N1), function(x)  { C=sum(N1[x,] %in% N2[x,]);  C / (2*k-C) }   ) )
    return(1-Jindex)
}
#e.g local.dist(  M1[study_annot=="facs",], fMNNcor[study_annot=="facs",], cosnorm=TRUE )



#### B2. Structure preservation. Global:
### Input are two K x d matrices M1, M2 with K samples and d features.
### n is the sample size of K
### The function first calculates two  n x n distance matrices corresponding to the n samples from M1, M2
### Then the vectorized forms of the two (upper/lower triangular) distance matrices are compared
### using one of: 1. KStest D statistic 2. 1-cosine similarity 3. 1-pearson's rho or 4. 1-Spearman's rho 
### to derive a single distance index.
### Ranges from 0-1. Lower is better.
library(coop)
library(scran)
global.dist <- function (m1, m2, metric=c("KStest","cosine","rho","Sp.rho"), n=1000, seed=1, cosnorm=FALSE ) {
    
    if (!identical(dim(m1), dim(m2))) {
        stop("The two input matrices must have the same dimensions")
    }
    
    
    if (seed){
        set.seed(seed)
    }
    
    s <- sample(1:nrow(m1), min(n, nrow(m1)) )
    m1  <- m1[s,]
    m2  <- m2[s,] 
    
    # Cosine normalization of feature matrices:
    if (cosnorm) {
        m1 <-   t(scran::cosineNorm( t(m1) ))
        m2 <-   t(scran::cosineNorm( t(m2) ))
    }
    
    #Compute distance matrices and vectorize:
    m1 <- as.matrix(dist(m1))
    m2 <- as.matrix(dist(m2))
    v1 <- as.vector(m1[upper.tri(m1)])
    v2 <- as.vector(m2[upper.tri(m2)])
    
    if (metric=="cosine"){
        D <- 1 - coop::cosine(v1 ,v2)
    }
    else if (metric=="rho"){
        D <- 1 - coop::pcor(v1, v2) 
    }
    else if (metric=="Sp.rho"){
        v1 <- rank(v1)
        v2 <- rank(v2)
        D <- 1 - coop::pcor(v1, v2)
    }
    else if (metric=="KStest"){
        D <- ks.test(v1 - mean(v1), v2 - mean(v2) )
        D <- D$statistic
    }
    return(D)
}

#e.g global.dist(  M1[study_annot=="facs",], fMNNcor[study_annot=="facs",], cosnorm=TRUE, n=2000 )



