#' Plot cells in the ebedded space
#'
#' Plot cells in the bidimensional space and color it according to a specific parameter.
#' 
#' @param data list; GFICF object
#' @param colorBy characters; Color cells according to a column contained in data$embedded data frame. Default is NULL.
#' @param pointSize integer; Size of the points in the plot. Default is 0.5.
#' @return The updated gficf object.
#' @export
#' @import ggrepel
#' @import ggplot2
#' 
#' @export
plotCells = function(data,colorBy=NULL,pointSize=.5)
{
  if (is.null(colorBy)) {return(ggplot(data = data$embedded) + geom_point(aes(x=X,y=Y),size=pointSize,color="blue") + theme_bw())}
  
  if (!colorBy%in%colnames(data$embedded)) {stop("colorBy parameter not found")}
  df = data$embedded
  u = unique(df[,colorBy])
  c = NULL
  for (i in u)
  {
    d = as.matrix(dist(df[df[,colorBy]%in%i,c(1,2)]))
    d = apply(d, 1, sum)
    ix = which.min(d)
    c = rbind(c,data.frame(cluster=i,xx=df[names(d[ix]),"X"],yy=df[names(d[ix]),"Y"],stringsAsFactors = F))
  }
  
  tmp = df[,colorBy]
  
  ggplot(data = data$embedded) + geom_point(aes(x=X,y=Y,color=tmp),size=pointSize) + theme_bw() + geom_text_repel(data = c,aes(x=xx,y=yy,label=cluster),min.segment.length = 0) + geom_point(data = c,aes(x=xx,y=yy),size=2) + theme(legend.position = "none")
}

#' Plot gene expression across cells
#'
#' Plot the expression of a group of genes across cells.
#' 
#' @param data list; GFICF object
#' @param genes characters; Id of genes to plot. It must correspond to the IDs on the rows of raw count matrix.
#' @param x Matrix; Custom normalized raw counts. If present will be used instead of the ones normalized by gficf. Default is NULL.
#' @return A list of plots.
#' @export
#' @import Matrix
#' @import ggplot2
#' 
#' @export
plotGenes = function(data,genes,x=NULL)
{
  if (is.null(data$embedded)) {stop("Please run reduction in the embedded space first!")}
  if (!is.null(x)) {data$rawCounts=x}
  if (is.null(data$rawCounts)) {stop("Raw or normalized counts absent.")}
  
  genes = genes[genes%in%rownames(data$rawCounts)]
  
  if (length(genes)==0) {stop("Genes are absent in the Expression matrix")}
  
  l  = vector(mode = "list",length = length(genes))
  names(l) = genes
  for (i in genes)
  {
    df = data$embedded
    df$expr = data$rawCounts[i,rownames(df)]
    df$expr <- (df$expr-min(df$expr))/(max(df$expr)-min(df$expr))
    df = df[order(df$expr,decreasing = F),]
    l[[i]] = ggplot(data = df,aes(x=X,y=Y,color=expr)) + geom_point(aes(shape=lev),size=.5,shape=19) + theme_bw() + theme_bw() + scale_color_gradient2(low = "gray",mid = "#2171b5",high = "#08306b",midpoint = .5) + ggtitle(i)
  }
  
  return(l)
}

