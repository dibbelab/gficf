#' Gene Frequency - Inverse Cell Frequency (GF-ICF)
#'
#' R implementation of the GF-ICF (https://www.frontiersin.org/articles/10.3389/fgene.2019.00734/abstract)
#' Thanks to 3’-end scRNA-seq approaches, we can now have an accurate estimation of gene expression without having to account for gene length,
#' thus the number of transcripts (i.e. UMI) associated to each gene, strictly reflects the frequency of a gene in a cell, exactly like a word in a document.
#' GFICF (Gene Frequency - Inverce Cell Frequency) is analugous of TF-IDF scoring method as defined for tex dada. With GFICF we consider a cell to be analogous
#' to a document, genes analogous to words and gene counts to be analogous of the word’s occurrence in a document.
#' 
#' @param M Matrix; UMI cell count matrix
#' @param cell_proportion_max integer; Remove genes present in more then to the specifided proportion (0,1). Default 1.
#' @param cell_proportion_min integer; Remove genes present in less then or equal to the specifided proportion (0,1). Default is 0.05 (i.e. 5 percent).
#' @param storeRaw logical; Store UMI counts.
#' @param  normalize logical; Rescale UMI counts before applay GFICF. Recaling is done using EdgeR normalization.
#' @param verbose boolean; Increase verbosity.
#' @return The updated gficf object.
#' @export
gficf = function(M,cell_proportion_max = 1,cell_proportion_min = 0.05,storeRaw=TRUE,normalize=TRUE,verbose=TRUE)
{
  data = list()
  M = normCounts(M,doc_proportion_max = cell_proportion_max,doc_proportion_min = cell_proportion_min,normalizeCounts=normalize,verbose=verbose)
  data$gficf = tf(M,verbose = verbose)
  if (storeRaw) {data$rawCounts=M;rm(M)}
  data$w = getIdfW(data$gficf,verbose = verbose)
  data$gficf = idf(data$gficf,data$w,verbose = verbose)
  data$gficf = t(l.norm(t(data$gficf),norm = "l2",verbose = verbose))
  gc()
  
  data$param <- list()
  data$param$cell_proportion_max = cell_proportion_max
  data$param$cell_proportion_min = cell_proportion_min
  data$param$normalized = normalize
  return(data)
}

#' @import Matrix
#' @importFrom edgeR DGEList calcNormFactors cpm
#' 
normCounts = function(M,doc_proportion_max = 1,doc_proportion_min = 0.01,normalizeCounts=FALSE,verbose=TRUE)
{
  ix = Matrix::rowSums(M!=0)
  M = M[ix>ncol(M)*doc_proportion_min & ix<=ncol(M)*doc_proportion_max,]
  
  if (normalizeCounts) 
  {
    tsmessage("Normalize counts..",verbose = verbose)
    M <- Matrix::Matrix(cpm(calcNormFactors(DGEList(counts=M),normalized.lib.sizes = T)),sparse = T) 
  } 
  
  return(M)
}


#' @import Matrix
#' 
tf = function(M,verbose)
{

  tsmessage("Apply GF transformation..",verbose = verbose)
  M =t(t(M) / Matrix::colSums(M))
  
  return(M)
}

#' @import Matrix
#' 
idf = function(M,w,verbose)
{
  tsmessage("Applay ICF..",verbose = verbose)
  M = M[rownames(M) %in% names(w),]
  if(nrow(M)<length(w))
  {
    g = names(w)[!names(w)%in%rownames(M)]
    tmp = Matrix::Matrix(data = 0,nrow = length(g),ncol = ncol(M))
    rownames(tmp) = g
    colnames(tmp) = colnames(M)
    M = rbind(M,tmp)
  }
  M = M[names(w),]
  M = M * w
  return(M)
}

#' @import Matrix
#' 
getIdfW = function(M,type="classic",verbose)
{
  tsmessage("Compute ICF weigth..",verbose = verbose)
  nt = Matrix::rowSums(M!=0)
  if (type == "classic") {w = log( (ncol(M)+1) / (nt+1) );rm(nt)}
  if (type == "prob") {w = log( (ncol(M) - nt) / nt );rm(nt)}
  if (type == "smooth") {w = log( 1 + ncol(M)/nt );rm(nt)}
  return(w)
}



l.norm = function (m, norm = c("l1", "l2"),verbose) 
{
  tsmessage(paste("Apply",norm),verbose = verbose)
  norm_vec = switch(norm, l1 = 1/rowSums(m), l2 = 1/sqrt(rowSums(m^2)))
  norm_vec[is.infinite(norm_vec)] = 0
  if (inherits(m, "sparseMatrix")) 
    Diagonal(x = norm_vec) %*% m
  else m * norm_vec
}

#' Save GFICF object.
#' 
#' Function to write a GFICF object to a file, preserving UMAP model.
#' 
#' @param data a GFICF object create by \code{\link{gficf}}.
#' @param file name of the file where the model is to be saved.
#' 
#' @examples 
#' # save
#' gficf_file <- tempfile("gficf_test")
#' saveGFICF(data, file = gficf_file)
#' 
#' # restore
#' data2 <- loadGFICF(file = gficf_file)
#' 
#' @export
saveGFICF <- function(data, file="~/Documents/gficf.test") 
{
  wd <- getwd()
  tryCatch({
    # create directory to store files in
    mod_dir <- tempfile(pattern = "dir")
    dir.create(mod_dir)
    
    # save gficf object
    gficf_dir <- file.path(mod_dir, "gficf")
    dir.create(gficf_dir)
    gficfl_tmpfname <- file.path(gficf_dir, "data")
    saveRDS(data, file = gficfl_tmpfname)
    
    # save uwot nn index if necessary
    if(data$reduction %in% c("umap","tumap"))
    {
      uwot_dir <- file.path(gficf_dir, "uwot_idx")
      dir.create(uwot_dir)
      nn_tmpfname <- file.path(uwot_dir,"nn1")
      data$uwot$nn_index$save(nn_tmpfname)
      data$uwot$nn_index$unload()
      data$uwot$nn_index$load(nn_tmpfname)
    }
    # archive the files under the temp dir into the single target file
    # change directory so the archive only contains one directory
    setwd(mod_dir)
    utils::tar(tarfile = file, files = "gficf/")
  },
  finally = {
    setwd(wd)
    if (file.exists(mod_dir)) {
      unlink(mod_dir, recursive = TRUE)
    }
  })
}

#' Restore GFICF object.
#' 
#' Function to read a GFICF object from a file saved with \code{\link{saveGFICF}}.
#' 
#' @param file name of the file where the object is stored.
#' 
#' @examples 
#' gficf_file <- tempfile("gficf_data")
#' 
#' # restore
#' data2 <- loadGFICF(file = gficf_file)
#' 
#' @export
loadGFICF <- function(file)
{
  model <- NULL
  
  tryCatch({
    # create directory to store files in
    mod_dir <- tempfile(pattern = "dir")
    dir.create(mod_dir)
    utils::untar(file, exdir = mod_dir)
    
    # load gficf object
    gficf_fname <- file.path(mod_dir, "gficf/data")
    if (!file.exists(gficf_fname)) {
      stop("Can't find gficf data in ", file)
    }
    data <- readRDS(file = gficf_fname)
    
    # load umap index if necessary
    if(data$reduction %in% c("umap","tumap"))
    {
      nn_fname <- file.path(mod_dir,"gficf/uwot_idx/nn1")
      if (!file.exists(nn_fname)) {
        stop("Can't find nearest neighbor index ", nn_fname, " in ", file)
      }
      
      # can provide any value for ndim as we get a new value when we load
      ann <- uwot:::create_ann(names(data$uwot$metric)[1], ndim = 1)
      ann$load(nn_fname)
      data$uwot$nn_index <- ann
    }
  },
  finally = {
    if (file.exists(mod_dir)) {
      unlink(mod_dir, recursive = TRUE)
    }
  })
  
  return(data)
}
