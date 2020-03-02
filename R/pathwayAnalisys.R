#' @importFrom  biomaRt useMart getBM getLDS
gmtPathways <- function(gmt.file,convertToEns,convertHu2Mm,verbose)
{
  pathwayLines <- strsplit(readLines(gmt.file), "\t")
  pathways <- lapply(pathwayLines, tail, -2)
  names(pathways) <- sapply(pathwayLines, head, 1)
  
  if (convertToEns & !convertHu2Mm)
  {
    tsmessage(".. Start converting human symbols to human ensamble id",verbose = verbose)
    g = as.character(unique(unlist(pathways)))
    ensembl = biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
    ens.map = biomaRt::getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),filters = 'hgnc_symbol',values = g,mart = ensembl,verbose = F)
    pathways = lapply(pathways, function(x,y=ens.map){r=y$ensembl_gene_id[y$hgnc_symbol%in%x];r=r[!is.na(r)];return(unique(r))})
    tsmessage("Done!",verbose = verbose)
  }
  
  if (convertToEns & convertHu2Mm )
  {
    tsmessage(".. Start converting human symbols to mouse ensamble id",verbose = verbose)
    g = as.character(unique(unlist(pathways)))
    human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    map.gene = biomaRt::getLDS(attributes = "hgnc_symbol",filters = "hgnc_symbol", values = g,mart = human,attributesL = "ensembl_gene_id", martL = mouse,verbose = F)
    pathways = lapply(pathways, function(x,y=map.gene) {r = unique(y$Gene.stable.ID[y$HGNC.symbol%in%x]);return(r[!is.na(r)])})
    tsmessage("Done!",verbose = verbose)
  }
  
  if (!convertToEns & convertHu2Mm )
  {
    tsmessage(".. Start converting human symbols to mouse symbols",verbose = verbose)
    g = as.character(unique(unlist(pathways)))
    human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    map.gene = biomaRt::getLDS(attributes = "hgnc_symbol",filters = "hgnc_symbol", values = g,mart = human,attributesL = "mgi_symbol", martL = mouse,verbose = F)
    pathways = lapply(pathways, function(x,y=map.gene) {r = unique(y$MGI.symbol[y$HGNC.symbol%in%x]);return(r[!is.na(r)])})
    tsmessage("Done!",verbose = verbose)
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
#' @param nt numeric; Number of cpu to use for the GSEA
#' @param minSize numeric; Minimal size of a gene set to test (default 15). All pathways below the threshold are excluded.
#' @param maxSize numeric; Maximal size of a gene set to test (default Inf). All pathways above the threshold are excluded.
#' @param verbose boolean; Show the progress bar.
#' @param seed integer; Seed to use for random number generation.
#' @return The updated gficf object.
#' @importFrom fgsea fgsea
#' @import fastmatch
#' @export
runGSEA <- function(data,gmt.file,nsim=1000,convertToEns=T,convertHu2Mm=F,nt=2,minSize=15,maxSize=Inf,verbose=TRUE,seed=180582)
{
  set.seed(seed)
  
  if (is.null(data$cluster.gene.rnk)) {stop("Please run clustcell function first")}
  
  data$gsea = list()
  data$gsea$pathways = gmtPathways(gmt.file,convertToEns,convertHu2Mm,verbose)
  data$gsea$es = Matrix::Matrix(data = 0,nrow = length(data$gsea$pathways),ncol = ncol(data$cluster.gene.rnk))
  data$gsea$nes = Matrix::Matrix(data = 0,nrow = length(data$gsea$pathways),ncol = ncol(data$cluster.gene.rnk))
  data$gsea$pval = Matrix::Matrix(data = 0,nrow = length(data$gsea$pathways),ncol = ncol(data$cluster.gene.rnk))
  data$gsea$fdr = Matrix::Matrix(data = 0,nrow = length(data$gsea$pathways),ncol = ncol(data$cluster.gene.rnk))
  
  rownames(data$gsea$es) = rownames(data$gsea$nes) = rownames(data$gsea$pval) = rownames(data$gsea$fdr) = names(data$gsea$pathways)
  colnames(data$gsea$es) = colnames(data$gsea$nes) = colnames(data$gsea$pval) = colnames(data$gsea$fdr) = colnames(data$cluster.gene.rnk)
  
  pb <- uwot:::Progress$new(max = ncol(data$cluster.gene.rnk),display = verbose)
  for (i in 1:ncol(data$cluster.gene.rnk))
  {
    df = as.data.frame(fgsea::fgsea(pathways = data$gsea$pathways,stats = data$cluster.gene.rnk[,i],nperm = nsim,gseaParam = 0,nproc = nt,minSize = minSize,maxSize = maxSize))[,1:7]
    data$gsea$es[df$pathway,i] = df$ES
    data$gsea$nes[df$pathway,i] = df$NES
    data$gsea$pval[df$pathway,i] = df$pval
    data$gsea$fdr[df$pathway,i] = df$padj
    pb$increment()
  }
  
  data$gsea$stat = df[,c("pathway","size")]
  return(data)
}

  