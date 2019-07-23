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