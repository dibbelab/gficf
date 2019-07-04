# gficf

An R implementation of the 
[Gene Frequency - Inverse Cell Frequency](https://link.to.the.paper.com) 
method for single cell data normalization (Gambardella et al. 2019), that also 
supports Phenograph Louvian Clustering using amazing annoy library from [uwot](https://github.com/jlmelville/uwot).
The package also support data reduction with either Principal Component Analisys (PCA) or
Latent Semantic Anlisys (LSA) before to apply t-SNE or UMAP for single cell data visualization.

## News

*July 3 2019*. First commit and documentation draft. 


## Installing

### From github

`gficf` makes use of annoy library in `uwot`. So you may have to carry out
a few extra steps before being able to build this package like for `uwot` installation:

**Windows**: install 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) and ensure 
`C:\Rtools\bin` is on your path.

**Mac OS X**: using a custom `~/.R/Makevars` 
[may cause linking errors](https://github.com/jlmelville/uwot/issues/1).
This sort of thing is a potential problem on all platforms but seems to bite
Mac owners more.
[The R for Mac OS X FAQ](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html#Installation-of-source-packages)
may be helpful here to work out what you can get away with. To be on the safe
side, I would advise building `uwot` without a custom `Makevars`.

```R
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
devtools::install_github("dibbelab/gficf")
```

## Example
Download Tabula Muris dataset from [HERE](https://drive.google.com/open?id=1yX8IQ7DiWG8PCmYieFFS7vj53Hf1OfT2) and annotation fro [HERE](https://drive.google.com/open?id=10ixOOsqZqf6GgwQP1okwoe_TMP_ZTzn5).

```R
library(gficf)

# See function man page for help
?gficf

# Common pipeline to use that goes from normalization to clustering

# Step 1: Nomrmalize data with gficf
data = gficf::gficf(M = readRDS("path/to/TabulaMuris.10x.mouse.RAW.rds"),cell_proportion_max = 1,cell_proportion_min = .05,storeRaw = F,normalize = F)

# Step 2: Reduce data with Latent Semantic Anlysis before to apply t-SNE or UMAP
data = gficf::runLSA(data = data,dim = 50)

# Alternative Step 2: Reduce data with Principal Component Analysis before to apply t-SNE or UMAP
# data = gficf::runPCA(data = data,dim = 50)

# Step 3: Applay t-SNE on reduced data and plot cells
data = gficf::runReduction(data = data,reduction = "tsne",seed = 0,nt=4)
gficf::plotCells(data = data)

# Alternative Step 3: Applay UMAP on reduced data
# data = gficf::runReduction(data = data,reduction = "umap",seed = 0,nt=4)
# gficf::plotCells(data = data)

# Step 4: Cell clustering using Phenograph approach
data = gficf::clustcells(data = data,from.embedded = F,dist.method = "manhattan",nt = 4,k = 50,community.algo = "louvian",seed = 0)

# Step 5: Visualize cells by identified clusters
gficf::plotCells(data = data,colorBy="cluster") + xlab("t-SNE1") + ylab("t-SNE2") + ggtitle("Cells colored by Clusters") 

```

![tabula_clusters.png](img/tabula_clusters.png) 

```R
# Additional steps: add annotation to cells and plot it.
info = readRDS("/path/to/TabulaMuris.10x.mouse.annotation.rds")
data$embedded$tissue = info$tissue
data$embedded$subtissue = info$subtissue
data$embedded$cell_ontology_class = info$cell_ontology_class
gficf::plotCells(data = data,colorBy="cell_ontology_class")

```
![tabula_annotated.png](img/tabula_annotated.png) 


