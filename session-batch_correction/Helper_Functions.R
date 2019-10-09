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





