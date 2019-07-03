#' Latent Semantic Analysis (LSA) 
#'
#' Reduce dimensionality of the single cell dataset using Latent Semantic Analysis (LSA)
#' 
#' @param data list; GFICF object
#' @param dim integer; Number of dimension which to reduce the dataset.
#' @param rescale logical; Rescale gficf scores before applying reduction (deprecated).
#' @param centre logical; Centre gficf scores before applying reduction (increase separation).
#' @param randomized logical; Use randomized (faster) version for matrix decomposition (default is TRUE).
#' @return The updated gficf object.
#' @importFrom RSpectra svds
#' @importFrom rsvd rsvd
#' 
#' @export
runLSA = function(data,dim=NULL,rescale=F,centre=F,randomized=T)
{
  if (is.null(dim))
  {
    if (is.null(data$dimPCA)) {stop("Specify the number of dims or run computePCADim first")} else {dim=data$dimPCA}
  } else {
    data$dimPCA = dim
  }
  
  x = t(data$gficf)
  if (rescale)
  {
    bc_tot <- Matrix::rowSums(x)
    median_tot <- stats::median(bc_tot)
    x <- sweep(x, 1, median_tot/bc_tot, '*')
  }
  
  if (centre)
  {
    x <- sweep(x, 2, Matrix::colMeans(x), '-')
    x <- sweep(x, 2, base::apply(x, 2, sd), '/')
  }
  
  if (randomized) {ppk<- rsvd::rsvd(x,k=dim)} else {ppk<- RSpectra::svds(x,k=dim)}
  data$pca <- ppk$u %*% base::diag(x = ppk$d)

  return(data)
}

#' Principal Component Analysis (PCA) 
#'
#' Reduce dimensionality of the single cell dataset using Principal Component Analysis (PCA)
#' 
#' @param data list; GFICF object
#' @param dim integer; Number of dimension which to reduce the dataset.
#' @param rescale logical; Rescale gficf scores before applying reduction (deprecated).
#' @param centre logical; Centre gficf scores before applying reduction (increase separation).
#' @param randomized logical; Use randomized (faster) version for matrix decomposition (default is TRUE).
#' @return The updated gficf object.
#' @importFrom rsvd rpca
#' 
#' @export
runPCA = function(data,dim=NULL,rescale=F,centre=F,randomized=T)
{
  if (is.null(dim))
  {
    if (is.null(data$dimPCA)) {stop("Specify the number of dims or run computePCADim first")} else {dim=data$dimPCA}
  } else {
    data$dimPCA = dim
  }
  
  x = data$gficf
  if (rescale)
  {
    bc_tot <- Matrix::colSums(x)
    median_tot <- stats::median(bc_tot)
    x <- sweep(x, 2, median_tot/bc_tot, '*')
  }
  
  if (centre)
  {
    x <- sweep(x, 1, Matrix::rowMeans(x), '-')
    x <- sweep(x, 1, base::apply(x, 1, sd), '/')
  }
  
  data$pca = rsvd::rpca(x,k=dim,center=F,scale=F,rand=randomized)$rotation
  
  return(data)
}

#' Dimensionality reduction
#'
#' Run t-SNE or UMAP or t-UMAP dimensionality reduction on selected features from PCA or LSA.
#' See ?umap or ?Rtsne for additional parameter to use. 
#' 
#' @param data list; GFICF object
#' @param reduction characters; Reduction method to use. One of:
#' \itemize{
#'   \item \code{"tsne"}
#'   \item \code{"umap"}
#'   \item \code{"tumap"} (the default)
#' }
#' @param nt integer; Number of thread to use (default 2).
#' @param seed integer; Initial seed to use.
#' @return The updated gficf object.
#' @import uwot
#' @importFrom Rtsne Rtsne
#' 
#' @export
runReduction = function(data,reduction="tumap",nt=2,seed=18051982, ...)
{

  reduction = base::match.arg(arg = reduction,choices = c("umap","tumap","tsne"),several.ok = F)
  
  set.seed(seed)
  if(reduction=="tumap")
  {
    data$embedded = base::as.data.frame(uwot::tumap(X = data$pca,scale = F,n_threads = nt,verbose = T, ...))
  }
  
  if(reduction=="umap")
  {
    data$embedded = base::as.data.frame(uwot::umap(X = data$pca, scale = F,n_threads = nt,verbose = T, ...))
  }
  
  if(reduction=="tsne")
  {
    data$embedded = base::as.data.frame(Rtsne::Rtsne(X = data$pca,dims = 2, pca = F,verbose = T,max_iter=1000,num_threads=nt, ...)$Y)
  }
  
  rownames(data$embedded) = base::colnames(data$gficf)
  colnames(data$embedded) = base::c("X","Y")
  return(data)
}

#' Number of features to use 
#'
#' Compute the number of dimension to use for either PCA or LSA.
#' 
#' @param data list; GFICF object
#' @param randomized logical; Use randomized (faster) version for matrix decomposition (default is TRUE).
#' @param subsampling logical; Use only a subset of the data for the imputation of dimensions to use.
#' @param plot logical; Show eblow plot.
#' @importFrom RSpectra svds
#' @importFrom rsvd rsvd
#' 
#' @export
computePCADim = function(data,randomized=T,subsampling=F,plot=T)
{
  dim = min(50,ncol(data$gficf))
  
  if (subsampling)
  {
    x = data$gficf[,sample(x = 1:ncol(data$gficf),size = round(ncol(data$gficf)/100*5))]
    if (randomized) {ppk<- rsvd::rsvd(t(x),k=dim)} else {ppk<- RSpectra::svds(t(x),k=dim)}
    rm(x)
  } else {
    if (randomized) {ppk<- rsvd::rsvd(t(data$gficf),k=dim)} else {ppk<- RSpectra::svds(t(data$gficf),k=dim)}
  }
  
  explained.var = ppk$d^2 / sum(ppk$d^2)
  if(plot) {plot(explained.var,xlab="components",ylab="explained.var")}
  
  ratio_to_first_diff <- diff(ppk$d^2 / sum(ppk$d^2)) / diff(ppk$d^2 / sum(ppk$d^2))[1]
  #reduction_dim <- (which(ratio_to_first_diff < 0.1) + 1)[1]
  ix = which(cumsum(diff((which(ratio_to_first_diff < 0.1))) == 1)>1)[1]
  reduction_dim = which(ratio_to_first_diff < 0.1)[ix]
  
  cat("Number of estimated dimensions =",reduction_dim)
  data$dimPCA = reduction_dim
  return(data)
}
