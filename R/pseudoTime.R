#' Pseudotime with slingshot
#'
#' The function is just a wrapper for slingshot and tradeSeq packages. 
#' Can be used to assign a pseudotime to each cell. Identified clusters in the UMAP space will be used.
#' 
#' @param data list; GFICF object
#' @param hgv boolean; use only top High Variable Genes (FDR<.25)
#' @param parallel boolean; Parralle execution
#' @param verbose boolean; Increase verbosity
#' @param start.cl integer; (optional) character, indicates the cluster(s) from which lineages will be drawn.
#' @param end.cl integer; (optional) character, indicates the cluster(s) which will be forced leaf nodes in their trees.
#' @param ... Additional arguments to pass to fitGAM of tradeSeq.
#' @return The updated gficf object.
#' @import slingshot
#' @import tradeSeq
#' @importFrom  BiocParallel multicoreWorkers
#' @importFrom  BiocParallel MulticoreParam
#' @importFrom  BiocParallel SerialParam
#' 
#' @export
runPseudoTime = function(data,hgv=T,parallel=F,verbose=T, start.cl = NULL, end.cl = NULL, ...)
{
  require(tradeSeq)
  if(is.null(data$community)) {stop("Pleas first run clustcells function to identify clusters!")}
  if(is.null(data$rawCounts)) {stop("Count data not stored!")}
  
  if(hgv) {
    tsmessage("... Computing HGV",verbose = verbose)
    df=findVarGenes(data = data$rawCounts,verbose = verbose)
    tsmessage("... Done!",verbose = verbose)
    g = df$gene[df$FDR<.25] # use top 25% var genes
  } else{
    g = rownames(data$rawCounts)
  }
  
  if(parallel) {
    cores = BiocParallel::multicoreWorkers()-1
    BiocParallel::register(BiocParallel::MulticoreParam(workers=cores))
  } else {
    BiocParallel::register(BiocParallel::SerialParam())
  }
  data$pseudoT = list()
  data$pseudoT$cl = data$embedded$cluster
  names(data$pseudoT$cl) = data$embedded$cell.id
  tsmessage("... Getting Linages",verbose = verbose)
  data$pseudoT$lin <- slingshot::getLineages(data = as.matrix(data$embedded[,c("X","Y")]),clusterLabels = data$pseudoT$cl,start.clus = start.cl, end.clus = end.cl)
  data$pseudoT$crv <- slingshot::getCurves(data$pseudoT$lin)
  tsmessage("... Done!",verbose = verbose)
  tsmessage("... (1/3) Fit GAM with tradeSeq",verbose = verbose)
  data$pseudoT$sce <- tradeSeq::fitGAM(counts = as.matrix(data$rawCounts[g,]),sds = data$pseudoT$crv,offset=rep(0,ncol(data$gficf)),parallel=parallel,verbose=verbose, ...)
  
  tsmessage("... (2/3) Computing DE genes across pseudotime (associationTest)",verbose = verbose)
  data$pseudoT$assoRes <- as.data.frame(tradeSeq::associationTest(data$pseudoT$sce))
  data$pseudoT$assoRes$fdr = p.adjust(data$pseudoT$assoRes$pvalue,method = "fdr")
  data$pseudoT$assoRes$gene.id = rownames(data$pseudoT$assoRes)
  
  tsmessage("... (3/3) Computing DE genes across pseudotime (startVsEndTest)",verbose = verbose)
  data$pseudoT$startRes <- as.data.frame(tradeSeq::startVsEndTest(data$pseudoT$sce))
  data$pseudoT$startRes$fdr = p.adjust(data$pseudoT$startRes$pvalue,method = "fdr")
  data$pseudoT$startRes$gene.id = rownames(data$pseudoT$startRes)
  tsmessage("... Done!",verbose = verbose)
  
  return(data)
}

make.gprege.param = function(n,ntimes)
{
  gpregeOptions <- list()
  gpregeOptions$indexRange=1:n
  gpregeOptions$explore <- FALSE
  gpregeOptions$exhaustPlotRes <- 30
  # Exhaustive plot contour levels.
  gpregeOptions$exhaustPlotLevels <- 10
  # Exhaustive plot maximum lengthscale.
  gpregeOptions$exhaustPlotMaxWidth <- 100
  # Noisy ground truth labels: which genes are in the top 786 ranks of the TSNI ranking.
  gpregeOptions$labels <-rep("FALSE",n)
  # SCG optimisation: maximum number of iterations.
  gpregeOptions$iters <- 100
  # SCG optimisation: no messages.
  gpregeOptions$display <- FALSE
  gpregeOptions$interpolatedT = 
    # Matrix of different hyperparameter configurations as rows:
    # [inverse-lengthscale percent-signal-variance percent-noise-variance].
    gpregeOptions$inithypers <- matrix( c( 1/1000, 1e-3, 0.999,1/ntimes, 0.999, 1e-3,1/80, 2/3, 1/3), ncol=3, byrow=TRUE)
  return(gpregeOptions)
}

#' @import gptk
get.gpscore = function(y,gp.opt)
{
  npsets <- dim(gp.opt$inithypers)[1]
  options <- gptk::gpOptions()
  options$kern$comp <- list("rbf", "white")
  loghypers <- matrix(0, 3, npsets)
  
  if (sum(is.nan(y)) > (length(y)/2))
  {
    message("Majority of points in profile are NaN.\n")
  }
  options$isMissingData <- any(is.nan(y))
  options$isSpherical <- !any(is.nan(y))
  stdy = sd(c(y[!is.nan(y)]))
  inithypers = t(log(gp.opt$inithypers %*% diag(c(1, stdy^2, stdy^2))))
  models <- list()
  LMLs = vector(mode = "numeric",length = npsets)
  for (h in 1:npsets)
  {
    models[[h]] <- gptk::gpCreate(1, 1, X = matrix(1:length(y),ncol=1), y = matrix(y,ncol=1), options)
    models[[h]] <- gptk::gpExpandParam(models[[h]], inithypers[,h])
    if (h != 1)
    {
      models[[h]] <- gptk::gpOptimise(models[[h]], gp.opt$display,gp.opt$iters)
      loghypers[, h] <- gptk::gpExtractParam(models[[h]],only.values = FALSE)
    }
    LMLs[h] <- gptk::gpLogLikelihood(models[[h]])
  }
  return(LMLs)
}

#' @importFrom BiocParallel multicoreWorkers
#' @import parallel
gprege <- function(data, gpregeOptions,parallel=T)
{
  gpregeOutput <- list()
  data <- t(as.matrix(scale(t(data), scale = FALSE)))
  
  if(parallel)
  {
    cl = parallel::makeCluster(BiocParallel::multicoreWorkers()-1)
    LMLs = parallel::parApply(cl = cl,X = data,MARGIN = 1,get.gpscore,gp.opt=gpregeOptions)
    parallel::stopCluster(cl)
  } else {
    LMLs = apply(data, 1, get.gpscore,gp.opt=gpregeOptions)
  }
  
  LMLs = t(LMLs)
  gpregeOutput$rankingScores = apply(as.matrix(LMLs[, 2:ncol(LMLs)]),1, max) - LMLs[, 1]
  return(gpregeOutput)
}

add.gp.score = function(data,curve=1)
{
  df = data$embedded
  df$pt = data$pseudoT$sce$slingshot[,paste0("pseudotime.curve",curve)]
  df = subset(df,!is.na(pt))
  int = data$pseudoT$sce@metadata$tradeSeq$knots
  df$pt_int[df$pt==0] = paste0("]",round(int[1],2),",",round(int[2],2),"]")
  lev = NULL
  for (i in 2:length(int)) {
    df$pt_int[df$pt>int[i-1] & df$pt<=int[i]] = paste0("]",round(int[i-1],2),",",round(int[i],2),"]")
    lev = c(lev,paste0("]",round(int[i-1],2),",",round(int[i],2),"]"))
  }
  df$pt_int = factor(as.character(df$pt_int),levels = lev)
  pt.mean = base::sapply(lev, function(x,y=log(data$rawCounts[data$pseudoT$assoRes$gene.id,df$cell.id]+1),z=as.character(df$pt_int)) Matrix::rowMeans(y[,z%in%x]))
  
  gpOptions <- make.gprege.param(n = nrow(pt.mean),ntimes = length(lev))
  output=suppressWarnings(gprege(data=as.matrix(pt.mean),gpregeOptions=gpOptions,parallel = T))
  
  data$pseudoT$assoRes[,paste0("gp.score.curve",curve)] = NA
  data$pseudoT$assoRes[,paste0("gp.score.curve",curve)] = output$rankingScores[rownames(data$pseudoT$assoRes)]
  data$pseudoT$startRes[,paste0("gp.score.curve",curve)] = NA
  data$pseudoT$startRes[,paste0("gp.score.curve",curve)] = output$rankingScores[rownames(data$pseudoT$startRes)]
  return(data)
}

plot.gene = function(data,ens,curve=1)
{
  df = data$embedded
  df$pt = data$pseudoT$sce$slingshot[,paste0("pseudotime.curve",curve)]
  int = data$pseudoT$sce@metadata$tradeSeq$knots
  df$pt_int[df$pt==0] = paste0("]",round(int[1],2),",",round(int[2],2),"]")
  lev = NULL
  for (i in 2:length(int)) {
    df$pt_int[df$pt>int[i-1] & df$pt<=int[i]] = paste0("]",round(int[i-1],2),",",round(int[i],2),"]")
    lev = c(lev,paste0("]",round(int[i-1],2),",",round(int[i],2),"]"))
  }
  df$expr = log10(data$rawCounts[ens,]+1)
  df$pt_int = factor(as.character(df$pt_int),levels = lev)
  stat = ddply(df,"pt_int",summarise,mu=mean(expr))
  ggplot(data=stat,aes(x=as.numeric(pt_int),y=mu)) + geom_point(size=3) + geom_smooth(method = lm, se = T) + theme_bw() + xlab("pseudotime") + ylab("log10(UMI+1)")
}


