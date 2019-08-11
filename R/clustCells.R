#' PhenoGraph clustering
#' 
#' R implementation of the PhenoGraph algorithm
#' 
#' A custom R implementation of the PhenoGraph (http://www.cell.com/cell/abstract/S0092-8674(15)00637-6) algorithm, 
#' which is a clustering method designed for high-dimensional single-cell data analysis. It works by creating a graph ("network") representing 
#' phenotypic similarities between cells by calclating the Jaccard coefficient between nearest-neighbor sets, and then identifying communities 
#' using the well known Louvain method (https://sites.google.com/site/findcommunities/) in this graph. 
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
#'   \item \code{"louvian"} (the default, the original Louvian method)
#'   \item \code{"louvian 2"} (Louvian with modularity optimization from Seurat)
#'   \item \code{"louvian 3"} (Louvain algorithm with multilevel refinement from Seurat)
#'   \item \code{"walktrap"}
#'   \item \code{"fastgreedy"}
#' }
#' @param store.graph logical; Store produced phenograph in the gficf object
#' @param seed integer; Seed to use for replication.
#' @param verbose logical; Increase verbosity.
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities (used only for louvian 2 or 3 methods).
#' @param n.start Number of random starts (used only for louvian 2 or 3 methods).
#' @param n.iter Maximal number of iterations per random start (used only for louvian 2 or 3 methods).
#' @return the updated gficf object
#' @importFrom  igraph graph.data.frame simplify cluster_louvain walktrap.community fastgreedy.community membership as_adj
#' @import uwot
#' @import Matrix
#' @export
clustcells <- function(data,from.embedded=F,k=15,dist.method="manhattan",nt=2,community.algo="louvian",store.graph=T,seed=180582,verbose=TRUE, resolution = 0.8, n.start = 10, n.iter = 10)
{
  community.algo = base::match.arg(arg = community.algo,choices = c("louvian","louvian 2","louvian 3","walktrap","fastgreedy"),several.ok = F)
  
  if (is.null(data$embedded)) {stop("Run first runReduction function")}
  set.seed(seed)
  if(verbose) {message("Finding Neighboors..")}
  
  if (from.embedded)
  {
    if(is.null(data$embedded)) {stop("First run runReduction to embed your cells")}
    neigh = uwot:::find_nn(as.matrix(data$embedded[,c(1,2)]),k=k+1,include_self = T,n_threads = nt,verbose = TRUE,method = "annoy",metric=dist.method)$idx
  } else {
    if(is.null(data$pca)) {stop("First run runPCA or runLSA to reduce dimensionality")}
    neigh = uwot:::find_nn(data$pca$cells,k=k+1,include_self = T,n_threads = nt,verbose = verbose,method = "annoy",metric=dist.method)$idx
  }
  
  neigh = neigh[,-1]
  relations <- jaccard_coeff(neigh,verbose)
  relations <- relations[relations[,1]>0, ]
  relations <- as.data.frame(relations)
  colnames(relations)<- c("from","to","weight")
  g <- igraph::graph.data.frame(relations, directed=FALSE)
  
  if (community.algo=="louvian")
  {
    if(verbose) {message("Performing louvain...")}
    community <- igraph::cluster_louvain(g)
  }
  
  if (community.algo=="louvian 2")
  {
    if(verbose) {message("Performing louvain with modularity optimization...")}
    community <- RunModularityClustering(igraph::as_adjacency_matrix(g,attr = "weight",sparse = T),1,resolution,1,n.start,n.iter,seed,verbose)
  }
  
  if (community.algo=="louvian 3")
  {
    if(verbose) {message("Performing louvain with modularity optimization...")}
    community <- RunModularityClustering(igraph::as_adjacency_matrix(g,attr = "weight",sparse = T),1,resolution,2,n.start,n.iter,seed,verbose)
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
  
  if(community.algo %in% c("louvian 2","louvian 3")) {
    community = community + 1
    data$embedded$cluster = as.character(community)
  } else {
    data$embedded$cluster <- as.character(igraph::membership(community))
  }
  
  if (store.graph) {data$community=community;data$cell.graph=g} else {data$community=community}
  
  # get centroid of clusters
  if(verbose) {message("Computing Centroids...")}
  cluster.map = data$embedded$cluster
  u = base::unique(cluster.map)
  data$cluster.centroids = base::sapply(u, function(x,y=data$gficf,z=cluster.map) Matrix::rowSums(y[,z%in%x]))
  
  message(paste("Detected Clusters:",length(unique(data$embedded$cluster))))
  
  return(data)
}

# Runs the modularity optimizer (C++ function from seurat package https://github.com/satijalab/seurat)
#
# @param SNN SNN matrix to use as input for the clustering algorithms
# @param modularity Modularity function to use in clustering (1 = standard; 2 = alternative)
# @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities
# @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python module.
# @param n.start Number of random starts
# @param n.iter Maximal number of iterations per random start
# @param random.seed Seed of the random number generator
# @param print.output Whether or not to print output to the console
# @param temp.file.location Deprecated and no longer used
# @param edge.file.name Path to edge file to use
#
# @return clusters
#
#' @importFrom utils read.table write.table
#
RunModularityClustering <- function(SNN = matrix(), modularity = 1, resolution = 0.8, algorithm = 1, n.start = 10, n.iter = 10, random.seed = 0, print.output = TRUE, temp.file.location = NULL, edge.file.name = "") 
{
  clusters <- RunModularityClusteringCpp(SNN,modularity,resolution,algorithm,n.start,n.iter,random.seed,print.output,edge.file.name)
  return(clusters)
}
