# gficf

An R implementation of the 
Gene Frequency - Inverse Cell Frequency method for single cell data
normalization [(Gambardella et al. 2019)](https://www.frontiersin.org/articles/10.3389/fgene.2019.00734/abstract).
The package also includes [Phenograph](https://www.cell.com/cell/fulltext/S0092-8674(15)00637-6)
[Louvain method](https://sites.google.com/site/findcommunities/)
clustering using [RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy) library
from [uwot](https://github.com/jlmelville/uwot).
The package also include data reduction with either Principal Component Analisys (PCA) or
Latent Semantic Analisys (LSA) before to apply t-SNE or UMAP for single cell data visualization.   

## News
*July 26 2019* **Added new functionality:** Is now possible to embed new cells in an already existing embedded space. [See example how..](https://jeky82.github.io/gficf_example.html#how-to-embedd-new-cells-in-an-existing-space)

*July 15 2019* Addedd web page with examples [HERE](https://jeky82.github.io/gficf_example.html)

*July 12 2019*. Paper Accepted and now available [HERE](https://www.frontiersin.org/articles/10.3389/fgene.2019.00734/abstract).

*July 3 2019*. Version 0.1 with example on Tabula Muris.


## Installing

### From github

`gficf` makes use of `uwot` and `Rcpp`. So you may have to carry out
a few extra steps before being able to build this package:

**Windows**: First install 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) and ensure 
`C:\Rtools\bin` is on your path.

**Mac OS X**: using a custom `~/.R/Makevars` may cause errors.
This sort of thing is a potential problem on all platforms but seems to bite
Mac owners more.
[The R for Mac OS X FAQ](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html#Installation-of-source-packages)
may be helpful here to work out what you can get away with. To be on the safe
side, I would advise building `gficf` without a custom `Makevars`.

```R
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
devtools::install_github("dibbelab/gficf")
```

## How to use GFICF
Web page with examples is now [HERE](https://jeky82.github.io/gficf_example.html)


## Useful Information

Apart from the man pages in R you may be interested in the following readings:

* A [description of t-SNE](https://lvdmaaten.github.io/tsne/).

* A [description of UMAP](https://jlmelville.github.io/uwot/umap-for-tsne.html)
using algorithmic terminology similar to t-SNE, rather than the more topological
approach of the UMAP publication.

* Some [Examples](https://jlmelville.github.io/uwot/umap-examples.html) of the 
output of UMAP on some datasets, compared to t-SNE. 

* Some results of running 
[UMAP on the simple datasets](https://jlmelville.github.io/uwot/umap-simple.html) 
from [How to Use t-SNE Effectively](https://distill.pub/2016/misread-tsne/).


## Phenograph Implementation Details
In the package `gficf` the function `clustcells` implement the [Phenograph](https://www.cell.com/cell/fulltext/S0092-8674(15)00637-6) algorithm,
which is a clustering method designed for high-dimensional single-cell data analysis. It works by creating a graph ("network") representing phenotypic similarities between cells by calculating the Jaccard coefficient between nearest-neighbor sets, and then identifying communities using the well known [Louvain method](https://sites.google.com/site/findcommunities/) in this graph. 

In this particular implementation of Phenograph we use approximate nearest neighbors found using [RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy)
libraries present in the `uwot` package. The supported distance metrics for KNN (set by the `dist.method` parameter) are:

* Euclidean (default)
* Cosine
* Manhattan
* Hamming

Please note that the Hamming support is a lot slower than the
other metrics. It is not recomadded to use it if you have more than a few hundred
features, and even then expect it to take several minutes during the index 
building phase in situations where the Euclidean metric would take only a few
seconds.

After computation of Jaccard distances among cells, the Louvain community detection is instead performed using `igraph` implementation.
All supported communities detection algorithm (set by the `community.algo` parameter) are:

* Louvain (default)
* Louvian with modularity optimization (c++ function imported from Seurat)
* Louvain algorithm with multilevel refinement (c++ function imported from Seurat)
* Walktrap
* Fastgreedy
