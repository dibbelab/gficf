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
  if (is.null(colorBy)) { return(ggplot(data = data$embedded) + geom_point(aes(x=X,y=Y),size=pointSize,color="blue") + theme_bw())}
  
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