#' Classify New Embedded Cells 
#'
#' Classify new embedded cells using GF-ICF transformation and K-nn algorithm.
#' Existing cells are used as training set.
#' 
#' @param data list; GFICF object
#' @param classes chareachters; Classes of aready exsiting cells in the order of thay are in colnames(data$gficf).
#' @param rescale k; Number of K-nn to use for classification. Odd number less than 30 are prefered.
#' @param seed integer; Initial seed to use.
#' @return A dataframe containing cell id and predicted classes.
#' @importFrom class knn
#' 
#' @export
classify.cells = function(data,classes,k=7,seed=18051982)
{
  set.seed(seed)
  classes = factor(as.character(classes))
  res = class::knn(data$pca$cells,data$pca$pred,classes,k = k,prob = F)
  df = data.frame(cell.id=rownames(data$pca$pred),pred=as.character(res),stringsAsFactors = F)
  return(df)
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
#' @param verbose boolean; Icrease verbosity.
#' @return The updated gficf object.
#' @import Matrix
#' @import uwot
#' @importFrom Rtsne Rtsne
#' 
#' @export
embedNewCells = function(data,x,nt=2,seed=18051982, verbose=TRUE, ...)
{
  x = normCounts(x[rownames(x)%in% names(data$w),],doc_proportion_max = 2,doc_proportion_min = 0,normalizeCounts = data$param$normalized,verbose=verbose)
  x = tf(x,verbose=verbose)
  x = idf(x,w = data$w,verbose=verbose)
  x = t(l.norm(t(x),norm = "l2",verbose=verbose))
  x = as.matrix(scaleMatrix(t(x), data$pca$rescale,data$pca$centre)) %*% data$pca$genes
  gc()
  
  if(data$reduction%in%c("tumap","umap")) {
    df = as.data.frame(uwot::umap_transform(as.matrix(x),data$uwot,verbose = verbose))
    rownames(df) = rownames(x)
    colnames(df) = c("X","Y")
  }
  
  if(data$reduction=="tsne") {
    warning("Not Fully supported!! With t-SNE only PCA/LSA components are predicted while t-SNE is re-run again!")
    set.seed(seed)
    df = base::as.data.frame(Rtsne::Rtsne(X = as.matrix(rbind(data$pca$cells,x)),dims = 2, pca = F,verbose = verbose,max_iter=1000,num_threads=nt, ...)$Y)
    rownames(df) = c(rownames(data$pca$cells),rownames(x))
    colnames(df) = c("X","Y")
    data$embedded[1:nrow(data$pca$cells),c("X","Y")] = df[1:nrow(data$pca$cells),c("X","Y")]
    df = df[rownames(x),]
  }
  
  if(is.null(data$embedded$predicted)) {data$embedded$predicted = "NO"}
  
  if (ncol(data$embedded)>2) {
    df[,colnames(data$embedded)[3:ncol(data$embedded)]] = NA
    df$predicted = "YES"
  } else {
    df$predicted = "YES"
  }
  data$embedded = rbind(data$embedded,df)
  data$pca$pred = x
  data$embedded$predicted = factor(as.character(data$embedded$predicted),levels = c("NO","YES"))
  return(data)
}


