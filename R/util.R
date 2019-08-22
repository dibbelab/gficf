#' @import Matrix
#' 
scaleMatrix = function(x,rescale,centre)
{
  if (rescale)
  {
    message("Rescaling..")
    bc_tot <- Matrix::rowSums(x)
    median_tot <- stats::median(bc_tot)
    x <- base::sweep(x, 1, median_tot/bc_tot, '*')
    message("Rescaling Done!")
  }
  
  if (centre)
  {
    message("Centering data..")
    x <- base::sweep(x, 2, Matrix::colMeans(x), '-')
    x <- base::sweep(x, 2, base::apply(x, 2, sd), '/')
    message("Centering Done!")
  }
  return(x)
}


stime <- function() {
  format(Sys.time(), "%T")
}

# message with a time stamp
tsmessage <- function(..., domain = NULL, appendLF = TRUE, verbose = TRUE,time_stamp = TRUE) {
  if (verbose) {
    msg <- ""
    if (time_stamp) {
      msg <- paste0(stime(), " ")
    }
    message(msg, ..., domain = domain, appendLF = appendLF)
    utils::flush.console()
  }
}

#' Convert Ensamble IDs to Official Gene Symbols
#'
#' It uses biomart. If more the one gene is associated to the enamble, the first one retrived from
#' Biomart is used.
#' 
#' @param df data frame; Data frame containing the IDs to convert.
#' @param col characters; Name of column containing the ensamble ids.
#' @param organism characters; Organism of origin (i.e. human or mouse).
#' @return The updated data frame with a new column called symb.
#'
#' @importFrom  biomaRt useMart getBM
#' @import fastmatch
#' 
#' @export
ensToSymbol = function(df,col,organism,verbose=T)
{
  organism = base::match.arg(arg = organism,choices = c("human","mouse"),several.ok = F)
  
  if(organism %in% "human")
  {
    tsmessage(".. Start converting human ens to human symbol",verbose = verbose)
    g = unique(as.character(df[,col]))
    ensembl = biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
    ens.map = biomaRt::getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),filters = 'ensembl_gene_id',values = g,mart = ensembl,verbose = F)
    df$symb = NA
    df$symb = ens.map$hgnc_symbol[fastmatch::fmatch(df[,col],ens.map$ensembl_gene_id)]
    tsmessage("Done!",verbose = verbose)
  }
  
  if (organism %in% "mouse")
  {
    tsmessage(".. Start converting human symbols to mouse ensamble id",verbose = verbose)
    g = as.character(unique(unlist(pathways)))
    ensembl = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    ens.map = biomaRt::getBM(attributes=c('ensembl_gene_id','mgi_symbol'),filters = 'ensembl_gene_id',values = g,mart = ensembl,verbose = F)
    df$symb = NA
    df$symb = ens.map$mgi_symbol[fastmatch::fmatch(df[,col],ens.map$ensembl_gene_id)]
    tsmessage("Done!",verbose = verbose)
  }
  
  return(df)
}
