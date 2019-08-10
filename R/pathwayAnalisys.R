#' @importFrom  biomaRt useMart getBM getLDS
gmtPathways <- function(gmt.file,convertToEns,convertHu2Mm)
{
  pathwayLines <- strsplit(readLines(gmt.file), "\t")
  pathways <- lapply(pathwayLines, tail, -2)
  names(pathways) <- sapply(pathwayLines, head, 1)
  
  if (convertToEns & !convertHu2Mm)
  {
    message(".. Start converting human symbols to human ensamble id")
    g = as.character(unique(unlist(pathways)))
    ensembl = biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
    ens.map = biomaRt::getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),filters = 'hgnc_symbol',values = g,mart = ensembl)
    pathways = lapply(pathways, function(x,y=ens.map){r=y$ensembl_gene_id[y$hgnc_symbol%in%x];r=r[!is.na(r)];return(unique(r))})
    message("Done!")
  }
  
  if (convertToEns & convertHu2Mm )
  {
    message(".. Start converting human symbols to mouse ensamble id")
    g = as.character(unique(unlist(pathways)))
    human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    map.gene = biomaRt::getLDS(attributes = "hgnc_symbol",filters = "hgnc_symbol", values = g,mart = human,attributesL = "ensembl_gene_id", martL = mouse)
    pathways = lapply(pathways, function(x,y=map.gene) {r = unique(y$Gene.stable.ID[y$HGNC.symbol%in%x]);return(r[!is.na(r)])})
    message("Done!")
  }
  
  if (!convertToEns & convertHu2Mm )
  {
    message(".. Start converting human symbols to mouse symbols")
    g = as.character(unique(unlist(pathways)))
    human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    map.gene = biomaRt::getLDS(attributes = "hgnc_symbol",filters = "hgnc_symbol", values = g,mart = human,attributesL = "mgi_symbol", martL = mouse)
    pathways = lapply(pathways, function(x,y=map.gene) {r = unique(y$MGI.symbol[y$HGNC.symbol%in%x]);return(r[!is.na(r)])})
    message("Done!")
  }
  
  pathways = pathways[sapply(pathways, length)>0]  
  pathways
}

#' Gene Set Enrichement Analysi on GF-ICF
#'
#' Compute GSEA for each cluster across a set of input pathways.
#' 
#' @param data list; GFICF object
#' @param gmt.file characters; Path to gmt file from MSigDB
#' @param nsim integer; number of simulation used to compute ES significance.
#' @param convertToEns boolean: Convert gene sets from gene symbols to Ensable id.
#' @param convertHu2Mm boolean: Convert gene sets from human symbols to Mouse Ensable id.
#' @return The updated gficf object.
#' @importFrom fgsea fgsea
#' @importFrom BiocParallel bpparam
#' @import fastmatch
#' @export
runGSEA <- function(data,gmt.file,nsim=1000,convertToEns=T,convertHu2Mm=F)
{
  if (is.null(data$cluster.centroids)) {stop("Please run clustcell function first")}
  
  data$gsea = list()
  data$gsea$pathways = gmtPathways(gmt.file,convertToEns,convertHu2Mm)
  data$gsea$es = Matrix::Matrix(data = 0,nrow = length(data$gsea$pathways),ncol = ncol(data$cluster.centroids))
  data$gsea$nes = Matrix::Matrix(data = 0,nrow = length(data$gsea$pathways),ncol = ncol(data$cluster.centroids))
  data$gsea$pval = Matrix::Matrix(data = 0,nrow = length(data$gsea$pathways),ncol = ncol(data$cluster.centroids))
  data$gsea$fdr = Matrix::Matrix(data = 0,nrow = length(data$gsea$pathways),ncol = ncol(data$cluster.centroids))
  
  rownames(data$gsea$es) = rownames(data$gsea$nes) = rownames(data$gsea$pval) = rownames(data$gsea$fdr) = names(data$gsea$pathways)
  colnames(data$gsea$es) = colnames(data$gsea$nes) = colnames(data$gsea$pval) = colnames(data$gsea$fdr) = colnames(data$cluster.centroids)
  
  pb <- txtProgressBar(min = 0, max = ncol(data$cluster.centroids), style = 3)
  for (i in 1:ncol(data$cluster.centroids))
  {
    df = as.data.frame(fgsea::fgsea(pathways = data$gsea$pathways,stats = data$cluster.centroids[,i],nperm = nsim,gseaParam = 0,BPPARAM = BiocParallel::bpparam()))[,1:7]
    data$gsea$es[,i] = df$ES
    data$gsea$nes[,i] = df$NES
    data$gsea$pval[,i] = df$pval
    data$gsea$fdr[,i] = df$padj
    setTxtProgressBar(pb, i)
  }
  
  data$gsea$stat = df[,c("pathway","size")]
  return(data)
}

  