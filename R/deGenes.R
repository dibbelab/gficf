#' Find Marker Genes from cell clusters.
#'
#' Try to identify marker genes across clusters performing Mann-Whitney U test.
#' DE genes are identified the expression in each cluster versus the all the other.
#'
#' @param data list; GFICF object
#' @param nt integer; Number of thread to use (default 2).
#' @param hvg boolean; Use only High Variable Genes (default is TRUE).
#' @param verbose boolean; Icrease verbosity.
#' @return The updated gficf object.
#' @import Matrix
#' @importFrom RcppParallel setThreadOptions
#' 
#' @export
findClusterMarkers = function(data,nt=2,hvg=T,verbose=T)
{
  if (is.null(data$community)) {stop("Please identify cluster first! Run clustcells function.")}
  if (is.null(data$rawCounts)) {stop("No raw/normalized counts stored. You should have run gficf normalization with storeRaw = T")}
  
  RcppParallel::setThreadOptions(numThreads = nt)
  
  u = unique(data$embedded$cluster)
  
  # Normalize if were not
  cpms = normCounts(data$rawCounts,doc_proportion_max = 2,doc_proportion_min = 0,normalizeCounts = !data$param$normalized,verbose=verbose)
  
  if (hvg)
  {
    tsmessage("... Detecting HVGs")
    df = findVarGenes(cpms,fitMethod = "locfit",verbose = verbose)
    df = df[order(df$FDR),]
    ix = df$FDR<.1
    tsmessage(paste("... Detected",sum(ix),"significant HVGs with an FDR < 10%"))
    if(sum(ix)<=1000) {ix = 1:min(1000,nrow(df))} 
    cpms = cpms[df$gene[ix],]
    rm(df)
  }
  
  cpms = as.matrix(cpms)
  
  tsmessage("... Start identify marker genes")
  res = NULL
  #progress_for(n=0,tot = length(u),display = T)
  for (i in 1:length(u)) {
    cells.1 = which(data$embedded$cluster%in%u[i])
    cells.2 = which(!data$embedded$cluster%in%u[i])
    tmp = rcpp_parallel_WMU_test(matX = cpms[,cells.1],matY = cpms[,cells.2],printOutput = F)
    tmp = data.frame(ens=rownames(cpms),log2FC=tmp[,2],p.value=tmp[,1],fdr=p.adjust(tmp[,1],method = "fdr"),stringsAsFactors = F)
    tmp = subset(tmp,fdr<.05 & log2FC>0)
    tmp = tmp[order(tmp$fdr,decreasing = F),]
    tmp$cluster = u[i]
    res = rbind(res,tmp)
    #progress_for(n=i,tot = length(u),display = T)
  }
  
  res = res[order(res$log2FC,decreasing = T),]
  rownames(res) = NULL
  data$de.genes = res
  return(data)
}

#'  Find high variable genes following the approach
#'  proposed by Chen et al. in BMC Genomics (2016)
#'  Code adapted from https://github.com/hillas/scVEGs
#' @param data list; GFICF object
#' @param fitMethod charachter; Method to use to fit variance and mean expression relationship (loess or locfit).
#' @param verbose boolean; Increase verbosity.
#' @import Matrix
#' @importFrom locfit locfit locfit.robust lp
#' @importFrom MASS fitdistr
#' 
findVarGenes = function(data,fitMethod="locfit",verbose=T)
{
  m <- dim(data)[1]
  std <- apply(data, 1, stats::sd)
  avg <- Matrix::rowMeans(data)
  cv <- std / avg
  # over dispersion sigma  (var = u(1 + u * sigma^2))
  xdata <- (avg)
  ydata <- log10(cv)
  xdata <- xdata[is.na(ydata) != "TRUE"]
  ydata <- ydata[is.na(ydata) != "TRUE"]
  
  if (fitMethod=="loess"){
    fitLoc <- stats::loess(formula = ydata ~ log10(x = xdata),data = data.frame("ydata"=ydata,"xdata"=xdata,stringsAsFactors = F),span = 0.8)
  } else {
    fitLoc <- locfit::locfit.robust(ydata ~ locfit::lp(log10(xdata), nn = .2))
  }
  
  xSeq <- seq(min(log10(xdata)), max(log10(xdata)), 0.005)
  gapNum <- matrix(0, length(xSeq), 1)
  for(i in 1:length(xSeq)) {
    cdx <- which((log10(xdata) >= xSeq[i] - 0.05) & (log10(xdata) < xSeq[i] + 0.05))
    gapNum[i,1] <- length(cdx)
  }
  cdx <- which(gapNum > m*0.005)
  xSeq <- 10 ^ xSeq
  ySeq <- predict(fitLoc,log10(xSeq))
  yDiff <- diff(ySeq)
  ix <- which(yDiff > 0 & log10(xSeq[-1]) > 0)
  if(length(ix) == 0)
    ix <- length(ySeq) - 1
  xSeq_all <- 10^seq(min(log10(xdata)), max(log10(xdata)), 0.001)
  xSeq <- xSeq[cdx[1]:ix[1] + 1]
  ySeq <- ySeq[cdx[1]:ix[1] + 1]
  
  b <- 1
  a <- 0
  df <- data.frame(x=xSeq, y = ySeq)
  fit = stats::nls(y ~ 0.5 * log10(b / x + a), data = df, start=list(b = b,a = a), stats::nls.control(maxiter = 500), na.action =  'na.exclude')
  newdf <- data.frame(x = xSeq_all)
  ydataFit <- stats::predict(fit,newdata = newdf)
  
  # Calculate CV difference
  logX <- log10(xdata)
  
  logXseq <- log10(xSeq_all)
  cvDist <- matrix(0,length(xdata),1)
  
  #progress_for(n=0,tot = length(logX),display = verbose)

  for (i in 1:length(logX)) 
  {
    cx <- which(logXseq >= logX[i] - 0.2 & logXseq < logX[i] + 0.2)
    tmp <- sqrt((logXseq[cx] - logX[i])^2 + (ydataFit[cx] - ydata[i])^2)
    tx <- which.min(tmp)
    
    if(logXseq[cx[tx]] > logX[i]) {
      if(ydataFit[cx[tx]] > ydata[i]) {
        cvDist[i] <- -1*tmp[tx]
      } else {
        cvDist[i] <- tmp[tx]
      }
      cvDist[i] <- -1*tmp[tx]
    } else if (logXseq[cx[tx]] <= logX[i]) {
      if(ydataFit[cx[tx]] < ydata[i]) {
        cvDist[i] <- tmp[tx]
      } else {
        cvDist[i] <- -1*tmp[tx]
      }
    }
    
    #progress_for(n=i,tot = length(logX),display = verbose)
    
  }
  
  cvDist <- log2(10^cvDist)
  
  # use kernel density estimate to find the peak
  dor <- stats::density(cvDist, kernel = "gaussian")
  distMid <-dor$x[which.max(dor$y)]
  dist2 <- cvDist - distMid
  tmpDist <- c(dist2[dist2 <= 0], abs(dist2[dist2 < 0])) + distMid
  distFit <- MASS::fitdistr(tmpDist, "normal")
  
  res = data.frame(gene=rownames(data),
                  mean=avg,
                  "cv"=cv,
                  P=pnorm(cvDist, mean = distFit$estimate[1], sd = distFit$estimate[2], lower.tail = FALSE)
  )
  res$FDR <- stats::p.adjust(res$P, 'fdr')
  
  return(res)
}

