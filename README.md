# gficf package overview

An R implementation of the 
Gene Frequency - Inverse Cell Frequency method for single cell data
normalization [(Gambardella et al. 2019)](https://www.frontiersin.org/articles/10.3389/fgene.2019.00734/abstract).
The package also includes [Phenograph](https://www.cell.com/cell/fulltext/S0092-8674(15)00637-6)
[Louvain method](https://sites.google.com/site/findcommunities/)
clustering using [RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy) library
from [uwot](https://github.com/jlmelville/uwot) and a naive but fast parallel implementation
of Jaccard Coefficient estimation using [RcppParallel](https://cran.r-project.org/package=RcppParallel).
The package also include data reduction with either Principal Component Analisys (PCA) or
Latent Semantic Analisys (LSA) before to apply t-SNE or UMAP for single cell data visualization.   

**Examples & Functionality**:
* General and simple use case scenario is [HERE](https://jeky82.github.io/gficf_example.html)

* Embed new cells in an already existing embedded space. [See example how..](https://jeky82.github.io/gficf_example.html#how-to-embedd-new-cells-in-an-existing-space)

* Idetify active pathways in a group of cells. [See example how..](https://jeky82.github.io/gficf_example.html#how-to-perform-gsea-to-identify-active-pathways-in-each-cluster)

* Idetify marker genes across clusters. [See example how..](https://jeky82.github.io/gficf_example.html#find-marker-genes)

## News
*Aug. 22 2019* **New functionality:** Identify marker genes across clusters of cells. [See example how..](https://jeky82.github.io/gficf_example.html#find-marker-genes)

*Aug. 20 2019* RcppParallel Mann–Whitney U test [(Benchmarks against R implementation)](https://jeky82.github.io/2019/08/20/MannWhitney.html) 

*Aug. 13 2019* **New functionality:** Identify active pathways in a group of cells. [See example how..](https://jeky82.github.io/gficf_example.html#how-to-perform-gsea-to-identify-active-pathways-in-each-cluster)

*Aug. 12 2019* RcppParallel Jaccard estimation in Phenograph [(20X speed boost with 6 cores)](https://jeky82.github.io/2019/08/12/parallel_JC_benchmarks.html) 

*Jul. 26 2019* **New functionality:** Embed new cells in an already existing embedded space. [See example how..](https://jeky82.github.io/gficf_example.html#how-to-embedd-new-cells-in-an-existing-space)

*Jul. 12 2019*. Paper Accepted and now available [HERE](https://www.frontiersin.org/articles/10.3389/fgene.2019.00734/abstract).

*Jul. 03 2019*. Version 0.1 with example on Tabula Muris.


## Installing

### From github

`gficf` makes use of `Rcpp`, `RcppParallel` and `RcppGSL`. So you may have to carry out
a few extra steps before being able to build this package:

**Linux systems**: You need gsl dev library to successfully install RcppGSL library.
On Ubuntu/Debian systems this can be accomplished by runnuing the command `sudo apt-get install libgsl-dev` from the terminal.

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

After computation of Jaccard distances among cells (custom [RcppParallel](https://cran.r-project.org/package=RcppParallel) implementation), the Louvain community detection is instead performed using `igraph` or native `Seurat`implementation.
All supported communities detection algorithm (set by the `community.algo` parameter) are:

* Louvain classic (default)
* Louvian with modularity optimization (native c++ function imported from `Seurat`)
* Louvain algorithm with multilevel refinement (native c++ function imported from `Seurat`)
* Walktrap
* Fastgreedy
