#' Gene Frequency - Inverse Document Frequency (GFICF)
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
#' @return The updated gficf object.
#' @export
gficf = function(M,cell_proportion_max = 1,cell_proportion_min = 0.05,storeRaw=TRUE,normalize=FALSE)
{
  data = list()
  M = gficf:::normCounts(M,doc_proportion_max = cell_proportion_max,doc_proportion_min = cell_proportion_min,normalizeCounts=normalize)
  data$gficf = gficf:::tf(M)
  data$w = gficf:::getIdfW(data$gficf)
  data$gficf = gficf:::idf(data$gficf,data$w)
  data$gficf = t(gficf:::l.norm(t(data$gficf),norm = "l2"))
  if (storeRaw) {data$rawCounts=M}
  
  data$param <- list()
  data$param$cell_proportion_max = cell_proportion_max
  data$param$cell_proportion_min = cell_proportion_min
  data$param$normalized = normalize
  return(data)
}

#' @import Matrix
#' @importFrom edgeR DGEList calcNormFactors cpm
#' 
normCounts = function(M,doc_proportion_max = 1,doc_proportion_min = 0.01,normalizeCounts=FALSE)
{
  ix = Matrix::rowSums(M!=0)
  M = M[ix>ncol(M)*doc_proportion_min & ix<=ncol(M)*doc_proportion_max,]
  
  if (normalizeCounts) 
  {
    message("Normalize counts..")
    M <- Matrix::Matrix(cpm(calcNormFactors(DGEList(counts=M),normalized.lib.sizes = T)),sparse = T) 
  } 
  
  return(M)
}


#' @import Matrix
#' 
tf = function(M)
{

  message("Apply GF transformation..")
  M =t(t(M) / Matrix::colSums(M))
  
  return(M)
}

#' @import Matrix
#' 
idf = function(M,w)
{
  message("Applay ICF..")
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
getIdfW = function(M,type="classic")
{
  message("Compute ICF weigth..")
  nt = Matrix::rowSums(M!=0)
  if (type == "classic") {w = log( (ncol(M)+1) / (nt+1) );rm(nt)}
  if (type == "prob") {w = log( (ncol(M) - nt) / nt );rm(nt)}
  if (type == "smooth") {w = log( 1 + ncol(M)/nt );rm(nt)}
  return(w)
}



l.norm = function (m, norm = c("l1", "l2")) 
{
  message(paste("Apply",norm))
  norm_vec = switch(norm, l1 = 1/rowSums(m), l2 = 1/sqrt(rowSums(m^2)))
  norm_vec[is.infinite(norm_vec)] = 0
  if (inherits(m, "sparseMatrix")) 
    Diagonal(x = norm_vec) %*% m
  else m * norm_vec
}
