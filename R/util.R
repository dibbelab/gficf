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
