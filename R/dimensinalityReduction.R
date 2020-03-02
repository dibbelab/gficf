#' Latent Semantic Analysis (LSA) 
#'
#' Reduce dimensionality of the single cell dataset using Latent Semantic Analysis (LSA)
#' 
#' @param data list; GFICF object
#' @param dim integer; Number of dimension which to reduce the dataset.
#' @param var.scale logical; Rescale gficf scores for adjusted variance like in pagoda2 (highly experimental!).
#' @param centre logical; Centre gficf scores before applying reduction (increase separation).
#' @param randomized logical; Use randomized (faster) version for matrix decomposition (default is TRUE).
#' @param seed integer; Initial seed to use.
#' @param use.odgenes boolean; Use significant overdispersed genes respect to ICF values (highly experimental!).
#' @param n.odgenes integer; Number of overdispersed genes to use (highly experimental!). A good choise seems to be usually between 1000 and 3000.
#' @param plot.odgenes boolean; Show significant overdispersed genes respect to ICF values.
#' @return The updated gficf object.
#' @importFrom RSpectra svds
#' @importFrom rsvd rsvd
#' 
#' @export
runLSA = function(data,dim=NULL,var.scale=F,centre=F,randomized=T,seed=180582,use.odgenes=F,n.odgenes=NULL,plot.odgenes=F)
{
  if(use.odgenes & is.null(data$rawCounts)) {stop("Raw Counts absent! Please run gficf normalization with storeRaw = T")}
  
  if (is.null(dim))
  {
    if (is.null(data$dimPCA)) {stop("Specify the number of dims or run computePCADim first")} else {dim=data$dimPCA}
  } else {
    data$dimPCA = dim
  }
  
  set.seed(seed)
  
  data$pca = list()
  data$pca$cells = t(data$gficf)
  #data$pca$cells = scaleMatrix(data$pca$cells,rescale,centre)
  
  if(use.odgenes) {
    overD=suppressWarnings(findOverDispersed(data = data,alpha = 0.1,verbose = F,plot = plot.odgenes))
    odgenes <- rownames(overD[overD$lpa<log(0.1),])
    if(!is.null(n.odgenes)) {
      if(n.odgenes>length(odgenes)) {
        odgenes <- rownames(overD)[(order(overD$lp,decreasing=F)[1:min(nrow(data$gficf),n.odgenes)])]
      } else {
        odgenes <- odgenes[1:n.odgenes]
      }
    }
    data$pca$cells = data$pca$cells[,odgenes]
    data$pca$odgenes = overD
    tsmessage("... using ",length(odgenes)," OD genes",verbose = T)
  } 
  
  if(var.scale) {
    if(!use.odgenes) {overD=suppressWarnings(findOverDispersed(data = data,alpha = 0.1,verbose = F))}
    data$pca$cells@x <- data$pca$cells@x*rep(overD[colnames(data$pca$cells),'gsf'],diff(data$pca$cells@p))
    data$pca$odgenes = overD
  }
  
  if (randomized) {ppk<- rsvd::rsvd(data$pca$cells,k=dim)} else {ppk<- RSpectra::svds(data$pca$cells,k=dim)}
  data$pca$cells <- ppk$u %*% base::diag(x = ppk$d)
  data$pca$centre <- F
  data$pca$rescale <- var.scale
  data$pca$genes <- ppk$v
  if(use.odgenes) {rownames(data$pca$genes)=odgenes} else {rownames(data$pca$genes) = rownames(data$gficf)}
  rownames(data$pca$cells) = colnames(data$gficf)
  return(data)
}

#' Principal Component Analysis (PCA) 
#'
#' Reduce dimensionality of the single cell dataset using Principal Component Analysis (PCA)
#' 
#' @param data list; GFICF object
#' @param dim integer; Number of dimension which to reduce the dataset.
#' @param var.scale logical; Rescale gficf scores for adjusted variance like in pagoda2 (highly experimental!).
#' @param centre logical; Centre gficf scores before applying reduction (increase separation).
#' @param randomized logical; Use randomized (faster) version for matrix decomposition (default is TRUE).
#' @param seed integer; Initial seed to use.
#' @param use.odgenes boolean; Use significant overdispersed genes respect to ICF values (highly experimental!).
#' @param n.odgenes integer; Number of overdispersed genes to use (highly experimental!). A good choise seems to be usually between 1000 and 3000.
#' @param plot.odgenes boolean; Show significant overdispersed genes respect to ICF values.
#' @return The updated gficf object. 
#' @return The updated gficf object.
#' @importFrom rsvd rpca
#' 
#' @export
runPCA = function(data,dim=NULL,var.scale=F,centre=F,randomized=T,seed=180582,use.odgenes=F,n.odgenes=NULL,plot.odgenes=F)
{
  
  if(use.odgenes & is.null(data$rawCounts)) {stop("Raw Counts absent! Please run gficf normalization with storeRaw = T")}
  
  if (is.null(dim))
  {
    if (is.null(data$dimPCA)) {stop("Specify the number of dims or run computePCADim first")} else {dim=data$dimPCA}
  } else {
    data$dimPCA = dim
  }
  
  set.seed(seed)
  
  data$pca = list()
  data$pca$cells = t(data$gficf)
  #data$pca$cells = scaleMatrix(data$pca$cells,F,F)
  
  if(use.odgenes) {
    overD=suppressWarnings(findOverDispersed(data = data,alpha = 0.1,verbose = F,plot = plot.odgenes))
    odgenes <- rownames(overD[overD$lpa<log(0.1),])
    if(!is.null(n.odgenes)) {
      if(n.odgenes>length(odgenes)) {
        odgenes <- rownames(overD)[(order(overD$lp,decreasing=F)[1:min(nrow(data$gficf),n.odgenes)])]
      } else {
        odgenes <- odgenes[1:n.odgenes]
      }
    }
    data$pca$cells = data$pca$cells[,odgenes]
    data$pca$odgenes = overD
    tsmessage("... using ",length(odgenes)," OD genes",verbose = T)
  }
  
  if(var.scale) {
    if(!use.odgenes) {overD=suppressWarnings(findOverDispersed(data = data,alpha = 0.1,verbose = F))}
    data$pca$cells@x <- data$pca$cells@x*rep(overD[colnames(data$pca$cells),'gsf'],diff(data$pca$cells@p))
    data$pca$odgenes = overD
  }
  
  x = rsvd::rpca(data$pca$cells,k=dim,center=centre,scale=F,rand=randomized)
  data$pca$cells = x$x
  data$pca$centre <- centre
  data$pca$rescale <- var.scale
  data$pca$genes <- x$rotation
  if(use.odgenes) {rownames(data$pca$genes)=odgenes} else {rownames(data$pca$genes) = rownames(data$gficf)}
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
#' @param verbose boolean; Icrease verbosity. 
#' @param ... Additional arguments to pass to Rtsne/umap/tumap call.
#' @return The updated gficf object.
#' @import uwot
#' @importFrom Rtsne Rtsne
#' 
#' @export
runReduction = function(data,reduction="tumap",nt=2,seed=18051982, ret_model_pred = T,verbose=T, ...)
{

  reduction = base::match.arg(arg = reduction,choices = c("umap","tumap","tsne"),several.ok = F)
  
  set.seed(seed)
  if (!is.null(data$pca))
  {
    if(reduction=="tumap"){
      data$uwot = uwot::tumap(X = data$pca$cells,scale = F,n_threads = nt,verbose = verbose,ret_model = ret_model_pred, ...)
      data$embedded = base::as.data.frame(data$uwot$embedding)
    }
    
    if(reduction=="umap"){
      data$uwot = uwot::umap(X = data$pca$cells, scale = F,n_threads = nt,verbose = verbose, ret_model = ret_model_pred, ...)
      data$embedded = base::as.data.frame(data$uwot$embedding)
    }
    
    if(reduction=="tsne"){
      data$uwot = NULL
      data$embedded = base::as.data.frame(Rtsne::Rtsne(X = data$pca$cells,dims = 2, pca = F,verbose = verbose,max_iter=1000,num_threads=nt, ...)$Y)
    }
  } else {
    message("Wrning: Reduction is applied directly on GF-ICF values.. can be slow if the dataset is big!")
    
    if(reduction=="tumap"){data$embedded = base::as.data.frame(uwot::tumap(X = as.matrix(t(data$gficf)),scale = F,n_threads = nt,verbose = verbose, ...))}
    
    if(reduction=="umap"){data$embedded = base::as.data.frame(uwot::umap(X = as.matrix(t(data$gficf)), scale = F,n_threads = nt,verbose = verbose, ...))}
    
    if(reduction=="tsne"){data$embedded = base::as.data.frame(Rtsne::Rtsne(X = as.matrix(t(data$gficf)), dims = 2, pca = F, verbose = verbose, max_iter=1000,num_threads=nt, ...)$Y)}
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

# find over dispersed genes respect to computed ICF
# ispired by pagoda2 function. Thanks to them.
#' @import mgcv
findOverDispersed=function(data,gam.k=5, alpha=5e-2, plot=FALSE, use.unadjusted.pvals=FALSE,do.par=T,max.adjusted.variance=1e3,min.adjusted.variance=1e-3,verbose=TRUE,min.gene.cells=0)
{
  rowSel <- NULL;
  
  tsmessage("calculating variance fit ...",verbose=verbose)
  df = data.frame(m=data$w,v=apply(data$rawCounts, 1, var),nobs=Matrix::rowSums(data$rawCounts!=0),stringsAsFactors = F)
  df$m = data$w
  
  
  # gene-relative normalizaton
  df$v <- log(df$v);
  rownames(df) <- rownames(data$gficf);
  vi <- which(is.finite(df$v) & df$nobs>=min.gene.cells);
  if(length(vi)<gam.k*1.5) { gam.k=1 };# too few genes
  if(gam.k<2) {
    tsmessage(" using lm ",verbose=verbose)
    m <- lm(v ~ m, data = df[vi,])
  } else {
    tsmessage(" using gam ",verbose=verbose)
    m <- mgcv::gam(as.formula(paste0('v ~ s(m, k = ',gam.k,')')), data = df[vi,])
  }
  df$res <- -Inf;  df$res[vi] <- resid(m,type='response')
  n.obs <- df$nobs;
  suppressWarnings(df$lp <- as.numeric(pf(exp(df$res),n.obs,n.obs,lower.tail=F,log.p=T)))
  df$lpa <- bh.adjust(df$lp,log=TRUE)
  n.cells <- ncol(data$gficf)
  df$qv <- as.numeric(qchisq(df$lp, n.cells-1, lower.tail = FALSE,log.p=TRUE)/n.cells)
  
  if(use.unadjusted.pvals) {
    ods <- which(df$lp<log(alpha))
  } else {
    ods <- which(df$lpa<log(alpha))
  }
  
  tsmessage(paste0(length(ods),'overdispersed genes ...',length(ods) ),verbose=verbose)
  
  df$gsf <- geneScaleFactors <- sqrt(pmax(min.adjusted.variance,pmin(max.adjusted.variance,df$qv))/exp(df$v));
  df$gsf[!is.finite(df$gsf)] <- 0;
  
  if(plot) {
    if(do.par) {
      par(mfrow=c(1,2), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
    }
    smoothScatter(df$m,df$v,main='',xlab='ICF value',ylab='log10[ variance ]')
    grid <- seq(min(df$m[vi]),max(df$m[vi]),length.out=1000)
    lines(grid,predict(m,newdata=data.frame(m=grid)),col="blue")
    if(length(ods)>0) {
      points(df$m[ods],df$v[ods],pch='.',col=2,cex=1)
    }
    smoothScatter(df$m[vi],df$qv[vi],xlab='ICF value',ylab='',main='adjusted')
    abline(h=1,lty=2,col=8)
    if(is.finite(max.adjusted.variance)) { abline(h=max.adjusted.variance,lty=2,col=1) }
    points(df$m[ods],df$qv[ods],col=2,pch='.')
  }
  tsmessage("done.\n",verbose=verbose)
  return(df)
}

# BH P-value adjustment with a log option
bh.adjust <- function(x, log = FALSE)
{
  nai <- which(!is.na(x))
  ox <- x
  x<-x[nai]
  id <- order(x, decreasing = FALSE)
  if(log) {
    q <- x[id] + log(length(x)/seq_along(x))
  } else {
    q <- x[id]*length(x)/seq_along(x)
  }
  a <- rev(cummin(rev(q)))[order(id)]
  ox[nai]<-a
  ox
}
