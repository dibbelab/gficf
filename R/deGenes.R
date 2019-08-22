#'  Find Marker Genes of Groups of Cells
#'
#' Try to identify marker genes across clusters performing Mann-Whitney U test
#' among each pair of clusters.
#'
#' @param data list; GFICF object
#' @param x Matrix; UMI counts matrix of cells to embedd.
#' @param nt integer; Number of thread to use (default 2).
#' @param hvg boolean; Use only High Variable Genes (default is TRUE).
#' @param verbose boolean; Icrease verbosity.
#' 
#' @return The updated gficf object.
#' @import Matrix
#' @importFrom RcppParallel setThreadOptions
#' 
#' @export
findClusterMarkers = function(data,nt=2,hvg=T,verbose=T)
{
  RcppParallel::setThreadOptions(numThreads = nt)
  
  u = unique(data$embedded$cluster)
  
  # Normalize if were not
  cpms = normCounts(data$rawCounts,doc_proportion_max = 2,doc_proportion_min = 0,normalizeCounts = !data$param$normalized,verbose=verbose)
  
  if (hvg)
  {
    tsmessage("... Detecting HVGs")
    df = findVarGenes(cpms,fitMethod = "locfit",verbose = verbose)
    tsmessage("... Done!")
    cpms = cpms[df$gene[df$FDR<.1],]
  }
  
  cpms = as.matrix(cpms)
  
  l = l2 = vector(mode = "list",length = length(u))
  names(l) = names(l2) = as.character(1:length(u))
  for (i in 1:length(l))
  {
    l[[i]] = Matrix::Matrix(data = 0,ncol = length(u)-1,nrow = nrow(cpms))
    colnames(l[[i]]) = names(l)[!names(l)%in%names(l)[i]]
    rownames(l[[i]]) = rownames(cpms)
    l2[[i]] = l[[i]]
  }
  
  tsmessage("... Identify DE genes among Cluster pairs")
  p = uwot:::Progress$new(max = round(.5 * length(u)^2),display = T)
  for (i in 1:length(u))
  {
    cells.1 = which(data$embedded$cluster%in%u[i])
    for (j in 1:i)
    {
      if(i!=j)
      {
        cells.2 = which(data$embedded$cluster%in%u[j])
        p_val <- gficf:::rcpp_parallel_WMU_test(matX = cpms[,cells.1],matY = cpms[,cells.2],F)
        l[[u[i]]][,u[j]] = p.adjust(p_val[,1],method = "fdr") * sign(p_val[,2])
        l[[u[j]]][,u[i]] = p.adjust(p_val[,1],method = "fdr") * sign(p_val[,2]) * -1
        
        l2[[u[i]]][,u[j]] = p_val[,2]
        l2[[u[j]]][,u[i]] = -1 * p_val[,2]
      }
      p$increment()
    }
  }
  
  tsmessage("... Identify Cluster markers")
  res = NULL
  for (i in 1:length(l))
  {
    M = l[[i]]
    M[abs(M)>=.5] = 0
    M = cbind(apply(M, 1, function(x) sum(x>0)),
              apply(M, 1, function(x) sum(x<0)),
              apply(M, 1, function(x) sum(x==0))
    )
    ix = (M[,1]==ncol(l[[i]])) #| M[,2]==ncol(l[[i]]))
    tmp = data.frame(ens = rownames(M)[ix],
                     fdr = apply(l[[i]][ix,],1,max),
                     fc  = Matrix::rowMeans(l2[[i]][ix,])
    )
    tmp$cluster = names(l)[i]
    res = rbind(res,tmp)
  }
  
  res = subset( res,!( ens %in% res$ens[duplicated(res$ens)] ) )
  rownames(res) = NULL
  data$de.genes = res
  return(res)
}

#'  Find high variable genes following the approach
#'  proposed by Chen et al. in BMC Genomics (2016)
#'  Code adapted from https://github.com/hillas/scVEGs
#'
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
  
  p = uwot:::Progress$new(max = length(logX),display = verbose)
  
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
    
    p$increment()
    
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

