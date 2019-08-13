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
runLSA = function(data,dim=NULL,rescale=F,centre=F,randomized=T,seed=180582)
{
  set.seed(seed)
  
  if (is.null(dim))
  {
    if (is.null(data$dimPCA)) {stop("Specify the number of dims or run computePCADim first")} else {dim=data$dimPCA}
  } else {
    data$dimPCA = dim
  }
  
  data$pca = list()
  data$pca$cells = t(data$gficf)
  data$pca$cells = scaleMatrix(data$pca$cells,rescale,centre)
  if (randomized) {ppk<- rsvd::rsvd(data$pca$cells,k=dim)} else {ppk<- RSpectra::svds(data$pca$cells,k=dim)}
  data$pca$cells <- ppk$u %*% base::diag(x = ppk$d)
  data$pca$centre <- centre
  data$pca$rescale <- rescale
  data$pca$genes <- ppk$v
  rownames(data$pca$genes) = rownames(data$gficf)
  rownames(data$pca$cells) = colnames(data$gficf)
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
runPCA = function(data,dim=NULL,rescale=F,centre=F,randomized=T,seed=180582)
{
  set.seed(seed)
  
  if (is.null(dim))
  {
    if (is.null(data$dimPCA)) {stop("Specify the number of dims or run computePCADim first")} else {dim=data$dimPCA}
  } else {
    data$dimPCA = dim
  }
  
  data$pca = list()
  data$pca$cells = t(data$gficf)
  data$pca$cells = scaleMatrix(data$pca$cells,rescale,centre)
  x = rsvd::rpca(data$pca$cells,k=dim,center=F,scale=F,rand=randomized)
  data$pca$cells = x$x
  data$pca$centre <- centre
  data$pca$rescale <- rescale
  data$pca$genes <- x$rotation
  rownames(data$pca$genes) = rownames(data$gficf)
  rownames(data$pca$cells) = colnames(data$gficf)
  colnames(data$pca$cells) = colnames(data$pca$genes) = paste("C",1:dim,sep = "")
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
#' @param ret_model_pred boolean; If true, the umap model is retained to be used for prediction. 
#' @param ... Additional arguments to pass to Rtsne/umap/tumap call.
#' @return The updated gficf object.
#' @import uwot
#' @importFrom Rtsne Rtsne
#' 
#' @export
runReduction = function(data,reduction="tumap",nt=2,seed=18051982, ret_model_pred = T, ...)
{

  reduction = base::match.arg(arg = reduction,choices = c("umap","tumap","tsne"),several.ok = F)
  
  set.seed(seed)
  if (!is.null(data$pca))
  {
    if(reduction=="tumap"){
      data$uwot = uwot::tumap(X = data$pca$cells,scale = F,n_threads = nt,verbose = T,ret_model = ret_model_pred, ...)
      data$embedded = base::as.data.frame(data$uwot$embedding)
    }
    
    if(reduction=="umap"){
      data$uwot = uwot::umap(X = data$pca$cells, scale = F,n_threads = nt,verbose = T, ret_model = ret_model_pred, ...)
      data$embedded = base::as.data.frame(data$uwot$embedding)
    }
    
    if(reduction=="tsne"){
      data$uwot = NULL
      data$embedded = base::as.data.frame(Rtsne::Rtsne(X = data$pca$cells,dims = 2, pca = F,verbose = T,max_iter=1000,num_threads=nt, ...)$Y)
    }
  } else {
    message("Wrning: Reduction is applied directly on GF-ICF values.. can be slow if the dataset is big")
    
    if(reduction=="tumap"){data$embedded = base::as.data.frame(uwot::tumap(X = as.matrix(t(data$gficf)),scale = F,n_threads = nt,verbose = T, ...))}
    
    if(reduction=="umap"){data$embedded = base::as.data.frame(uwot::umap(X = as.matrix(t(data$gficf)), scale = F,n_threads = nt,verbose = T, ...))}
    
    if(reduction=="tsne"){data$embedded = base::as.data.frame(Rtsne::Rtsne(X = as.matrix(t(data$gficf)),dims = 2, pca = F,verbose = T,max_iter=1000,num_threads=nt, ...)$Y)}
  }  
  rownames(data$embedded) = base::colnames(data$gficf)
  colnames(data$embedded) = base::c("X","Y")
  data$reduction = reduction
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

#' Embed new cells in an existing space 
#'
#' This function embed new cells in an already existing space. For now it supports only UMAP and t-UMAP. Briefly new cells are first normalized with GF-ICF method but using as ICF weigth estimated on the existing cells and than projected in the existing PCA/LSA space before to be embedded in the already existing UMAP space via umap_transform function. 
#' 
#' @param data list; GFICF object
#' @param x Matrix; UMI counts matrix of cells to embedd.
#' @param nt integer; Number of thread to use (default 2).
#' @param seed integer; Initial seed to use.
#' @param ... Additional arguments to pass to Rtsne or umap_transform call.
#' @return The updated gficf object.
#' @import Matrix
#' @import uwot
#' @importFrom Rtsne Rtsne
#' 
#' @export
embedNewCells = function(data,x,nt=2,seed=18051982, ...)
{
  x = normCounts(x[rownames(x)%in% names(data$w),],doc_proportion_max = 2,doc_proportion_min = 0,normalizeCounts = data$param$normalized)
  x = tf(x)
  x = idf(x,w = data$w)
  x = t(l.norm(t(x),norm = "l2"))
  pcapred = scaleMatrix(t(x), data$pca$rescale,data$pca$centre) %*% data$pca$genes
  rownames(pcapred) = colnames(x)
  colnames(pcapred) = colnames(data$pca$cells)
  
  if(data$reduction%in%c("tumap","umap")) {
    df = as.data.frame(uwot::umap_transform(as.matrix(pcapred),data$uwot,verbose = TRUE))
    rownames(df) = rownames(pcapred)
    colnames(df) = c("X","Y")
  }
  
  if(data$reduction=="tsne") {
    warning("Not Fully supported!! With t-SNE only PCA/LSA components are predicted while t-SNE is re-run again!")
    set.seed(seed)
    df = base::as.data.frame(Rtsne::Rtsne(X = as.matrix(rbind(data$pca$cells,pcapred)),dims = 2, pca = F,verbose = T,max_iter=1000,num_threads=nt, ...)$Y)
    rownames(df) = c(rownames(data$pca$cells),rownames(pcapred))
    colnames(df) = c("X","Y")
    data$embedded[1:nrow(data$pca$cells),c("X","Y")] = df[1:nrow(data$pca$cells),c("X","Y")]
    df = df[rownames(pcapred),]
  }
  
  if(is.null(data$embedded$predicted)) {data$embedded$predicted = "NO"}
  
  if (ncol(data$embedded)>2) {
    df[,colnames(data$embedded)[3:ncol(data$embedded)]] = NA
    df$predicted = "YES"
  } else {
    df$predicted = "YES"
  }
  data$embedded = rbind(data$embedded,df)
  data$pca$pred = pcapred
  data$embedded$predicted = factor(as.character(data$embedded$predicted),levels = c("NO","YES"))
  return(data)
}



