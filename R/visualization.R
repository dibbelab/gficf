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
#' @param log2Expr boolean; Relative expression of a gene is computed on rescaled in log2 expression (default TRUE).
#' @param x Matrix; Custom normalized raw counts. If present will be used instead of the ones normalized by gficf. Default is NULL.
#' @return A list of plots.
#' @import Matrix
#' @import ggplot2
#' 
#' @export
plotGenes = function(data,genes,log2Expr=T,x=NULL)
{
  if (is.null(data$embedded)) {stop("Please run reduction in the embedded space first!")}
  if (!is.null(x)) {data$rawCounts=x}
  if (is.null(data$rawCounts)) {stop("Raw or normalized counts absent.")}
  
  data$rawCounts = normCounts(data$rawCounts,doc_proportion_max = 2,
                              doc_proportion_min = 0,
                              normalizeCounts = !data$param$normalized & is.null(x),
                              verbose=T)
  
  genes = genes[genes%in%rownames(data$rawCounts)]
  
  if (length(genes)==0) {stop("Genes are absent in the Expression matrix")}
  
  l  = vector(mode = "list",length = length(genes))
  names(l) = genes
  for (i in genes)
  {
    df = data$embedded
    if(log2Expr)
    {
      df$expr = log2(data$rawCounts[i,rownames(df)]+1)
    } else {
      df$expr = data$rawCounts[i,rownames(df)]
    }
    
    df$expr = (df$expr-min(df$expr))/(max(df$expr)-min(df$expr))
    df = df[order(df$expr,decreasing = F),]
    l[[i]] = ggplot(data = df,aes(x=X,y=Y,color=expr)) + geom_point(aes(shape=lev),size=.5,shape=19) + theme_bw() + theme_bw() + scale_color_gradient2(low = "gray",mid = "#2171b5",high = "#08306b",midpoint = .5) + ggtitle(i)
    
  }
  
  return(l)
}

#' Plot the expression of a gene across group of cells.
#'
#' Plot the expression of a gene across group of cells with violion plot.
#' 
#' @param data list; GFICF object
#' @param genes characters; Id of genes to plot. It must correspond to the IDs on the rows of raw count matrix.
#' @param ncol numeric; Number of columns of the final plot (defaul is 3).
#' @param x Matrix; Custom normalized raw counts. If present will be used instead of the ones normalized by gficf. Default is NULL.
#' @return A list of plots.
#' 
#' @import fastmatch
#' @import Matrix
#' @import ggplot2
#' @importFrom reshape2 melt
#' 
#' @export
plotGeneViolin = function(data,gene,ncol=3,x=NULL)
{
  if (is.null(data$community)) {stop("Please run clustcells first!")}
  if (!is.null(x)) {data$rawCounts=x}
  if (is.null(data$rawCounts)) {stop("Raw or normalized counts absent.")}
  
  cpms = normCounts(data$rawCounts,doc_proportion_max = 2,
                    doc_proportion_min = 0,
                    normalizeCounts = !data$param$normalized & is.null(x),
                    verbose=T)
  
  df = reshape2::melt(as.matrix(cpms[gene,]))
  colnames(df) = c("ens","cell.id","value")
  
  if(!is.null(names(gene))) {
      ix = is.na(names(gene)) | names(gene) %in% "" | is.null(names(gene))
      if(sum(ix)>0) {names(gene)[ix] = gene[ix]}
      if(length(unique(names(gene))) == length(gene)) {
        df$ens = names(gene)[fastmatch::fmatch(df$ens,gene)]
      }
    }
  
  df$value = log2(df$value+1)
  df$cluster = data$embedded$cluster[match(df$cell.id,rownames(data$embedded))]
  df$cluster = factor(as.character(df$cluster),levels = as.character(1:length(unique(df$cluster))))
  p = ggplot2::ggplot(data = df,ggplot2::aes(x=cluster,y=value)) + 
      ggplot2::geom_violin(scale = "width") + 
      ggplot2::facet_wrap(~ens,scales = "free_y",ncol=ncol) + 
      ggplot2::ylab("log2(CPM+1)") + ggplot2::xlab("") + ggplot2::theme_bw()
  
  return(p)
}

#' Plot GSEA results
#'
#' Circle plot for gene set enrichement analysis results.
#' 
#' @param data list; GFICF object
#' @param fdr number; FDR threshold to select significant pathways to plot.
#' @return plot from ggplot2 package.
#' @import Matrix
#' @import ggplot2
#' @importFrom reshape2 melt
#' 
#' @export
plotGSEA = function(data,fdr=.05)
{
  if (is.null(data$gsea)) {stop("Please run runGSEA function first")}
  nes = data$gsea$nes
  nes[data$gsea$es<=0 | data$gsea$fdr>=fdr] = 0 
  nes = nes[Matrix::rowSums(nes)>0,]
  df = reshape2::melt(as.matrix(nes))
  colnames(df) = c("pathway","cluster","nes")
  ggplot(data = df,aes(x=pathway,y=cluster)) + geom_point(aes(size=nes)) + scale_size_continuous(range = c(0,7)) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_y_continuous(breaks = 1:max(df$cluster)) + xlab("") + ylab("Cluster name")
}

#' Plot GSEA results
#'
#' Circle plot for gene set enrichement analysis results.
#' 
#' @param data list; GFICF object
#' @param pathwayName characters; Name of the pathway to plot.
#' @param fdr number; FDR threshold to select significant pathways to plot.
#' @return plot from ggplot2 package.
#' @import Matrix
#' @import ggplot2
#' 
#' @export
plotPathway = function(data,pathwayName,fdr=.05)
{
  if (is.null(data$gsea)) {stop("Please run runGSEA function first")}
  nes = data$gsea$nes
  nes[data$gsea$es<=0 | data$gsea$fdr>=fdr] = 0 
  nes = nes[Matrix::rowSums(nes)>0,]
  nes = nes[pathwayName,]
  df = data$embedded
  df$NES = nes[match(df$cluster,names(nes))]
  ggplot(data = df,aes(x=X,y=Y)) + geom_point(aes(color=NES),shape=20) + theme_bw() + scale_color_gradient(low = "gray",high = "red")
}


