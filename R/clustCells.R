#' RphenoGraph clustering like
#' 
#' R implementation of the PhenoGraph algorithm
#' 
#' A custom R implementation of the PhenoGraph(http://www.cell.com/cell/abstract/S0092-8674(15)00637-6) algorithm, 
#' which is a clustering method designed for high-dimensional single-cell data analysis. It works by creating a graph ("network") representing 
#' phenotypic similarities between cells by calclating the Jaccard coefficient between nearest-neighbor sets, and then identifying communities 
#' using the well known Louvain method(https://sites.google.com/site/findcommunities/) in this graph. 
#' 
#' That version used PCA or LSA reduced meta-cells and multithreading annoy version for K-nn search (from uwot package).  
#' 
#' @param data list; Input data (gficf object)
#' @param from.embedded logical; Use embeddedd (UMAP or tSNA) space for clustering cells. Best results are usually obtained not using the embedded space.
#' @param k integer; number of nearest neighbours (default:15)
#' @param dist.method character; Dist to use for K-nn. Type of distance metric to use to find nearest neighbors. One of:
#' \itemize{
#'   \item \code{"euclidean"} (the default)
#'   \item \code{"cosine"}
#'   \item \code{"manhattan"}
#'   \item \code{"hamming"} (very slow)
#' }
#' @param  nt integer; Number of cpus to use for k-nn search
#' @param community.algo characthers; Community algorithm to use for clustering. Supported are:
#' \itemize{
#'   \item \code{"louvian"} (the default)
#'   \item \code{"walktrap"}
#'   \item \code{"fastgreedy"}
#' }
#' @param store.graph logical; Store produced phenograph in the gficf object
#' @param seed integer; Seed to use for replication.
#' @param verbose logical; Increase verbosity.
#' @return the updated gficf object
#' @importFrom  igraph graph.data.frame simplify cluster_louvain walktrap.community fastgreedy.community
#' @import uwot
#' @import Matrix
#' @export
clustcells <- function(data,from.embedded=F,k=15,dist.method="manhattan",nt=2,community.algo="louvian",store.graph=F,seed=180582,verbose=TRUE)
{
  set.seed(seed)
  if(verbose) {message("Finding Neighboors..")}
  
  if (from.embedded)
  {
    if(is.null(data$embedded)) {stop("First run runReduction to embed your cells")}
    neigh = uwot:::find_nn(as.matrix(data$embedded[,c(1,2)]),k=k,include_self = F,n_threads = nt,verbose = TRUE,method = "annoy",metric=dist.method)
  } else {
    if(is.null(data$pca)) {stop("First run runPCA or runLSA to reduce dimensionality")}
    neigh = uwot:::find_nn(data$pca,k=k,include_self = F,n_threads = nt,verbose = verbose,method = "annoy",metric=dist.method)
  }
  
  if(verbose) {message("Jaccard Coefficient..")}
  links <- jaccard_coeff(neigh$idx)
  links <- links[links[,1]>0, ]
  links <- links[links[,1] != links[,2],]
  relations <- as.data.frame(links)
  colnames(relations)<- c("from","to","weight")
  g <- igraph::graph.data.frame(relations, directed=FALSE)
  
  #sum multiple edges
  g = igraph::simplify(g, edge.attr.comb=list(weight="sum"),remove.loops = F,remove.multiple = T)
  
  if (community.algo=="louvian")
  {
    if(verbose) {message("Performing louvain...")}
    community <- igraph::cluster_louvain(g)
  }
  if (community.algo=="walktrap")
  {
    if(verbose) {message("Performing walktrap...")}
    community <- igraph::walktrap.community(g)
  } 
  if (community.algo=="fastgreedy") {
    if(verbose) {message("Performing fastgreedy...")}
    community <- igraph::fastgreedy.community(g)
  }
  
  data$embedded$cluster <- as.character(igraph::membership(community))
  
  if (store.graph) {data$community=community;data$cell.graph=g} else {data$community=community}
  
  # get centroid of clusters
  if(verbose) {message("Computing Centroids...")}
  cluster.map = as.character(membership(data$community))
  u = base::unique(cluster.map)
  data$cluster.centroids = base::sapply(u, function(x,y=data$gficf,z=cluster.map) Matrix::rowSums(y[,z%in%x]))
  
  message(paste("Detected Clusters:",length(unique(data$embedded$cluster))))
  
  return(data)
}

