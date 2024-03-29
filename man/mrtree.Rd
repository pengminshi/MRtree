% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrtree.R
\name{mrtree}
\alias{mrtree}
\alias{mrtree.matrix}
\alias{mrtree.data.frame}
\alias{mrtree.SingleCellExperiment}
\alias{mrtree.Seurat}
\title{Muti-resolution reconciled tree (MRtree)}
\usage{
mrtree(x, ...)

\method{mrtree}{matrix}(
  x,
  prefix = NULL,
  suffix = NULL,
  max.k = Inf,
  consensus = FALSE,
  sample.weighted = FALSE,
  augment.path = FALSE,
  verbose = FALSE,
  n.cores = 1
)

\method{mrtree}{data.frame}(x, prefix = NULL, suffix = NULL, ...)

\method{mrtree}{SingleCellExperiment}(x, prefix = "sc3_", suffix = "_clusters", ...)

\method{mrtree}{Seurat}(x, prefix = "RNA_snn_res.", suffix = NULL, ...)
}
\arguments{
\item{x}{Seurat object}

\item{...}{other parameters}

\item{prefix}{srting indicating columns containing clustering information}

\item{suffix}{string at the end of column names containing clustering information}

\item{max.k}{the maximum resolution (number of clusters) to consider in building the tree}

\item{consensus}{boolean, whether to perform consensus clustering within the clusterings with the same number of clusters.}

\item{sample.weighted}{boolean, whether to weight the samples by the size of clusters, where higher weight is given to rare clusters. False by default.}

\item{augment.path}{boolean, whether to augment one NA cluster in each layer to avoid shrinking the width of the tree with clustering alternate across layers. This will dramatically increase the running time. False by default.}

\item{verbose}{boolean, whether to show running messages. False by default.}

\item{n.cores}{integer, number of cores for parallel computation. 1 core by default.}
}
\value{
A list containing  \describe{
  \item{labelmat.mrtree}{The Reconciled tree saved as a label matrix, with duplicated layers omited.}
  \item{labelmat.recon}{The full reconciled tree save in a label matrix}
  \item{labelmat.flat}{The initial flat clustering (cluster tree) as input for MRtree algorithm}
  \item{resolutions}{The corresponding clustering resolution of the initial cluster tree}
  \item{paths}{The unique path in the resulted reconciled hierarchical cluster tree}
}
}
\description{
The method take as input the flat clusterings obtained across multiple resolution parameters. It reconciles the clusterings to produce a hierarchical cluster tree structure.

\code{mrtree} with label saved in a matrix as input

\code{mrtree} with labelmatrix save as a data frame as input

\code{mrtree} with SingleCellExperiment object as input

\code{mrtree} with Seurat object as input
}
\examples{
data("clust_example")
# matrix as input
out = mrtree(clust_example$clusterings)
# out$labelmat.mrtree
# data frame with given prefix and suffix as input
df = as.data.frame(clust_example$clusterings); colnames(df) = paste0('K_', 1:4)
df$other = 1:nrow(clust_example$clusterings)  # add an additional column
out = mrtree(df, prefix = 'K_', suffix=NULL)

cl = cbind(clust_example$clusterings, clust_example$clusterings)
# add some additional noise
for (i in 1:ncol(cl)){
    cl[sample(10),i] = sample(1:length(unique(cl[,i])), 10, replace = TRUE)
}
mrtree(cl, consensus=TRUE) # mrtree with within-resolution consensus clustering
mrtree(cl, augment.path=TRUE) # mrtree with augmentation of NA cluster
mrtree(cl, sample.weighted=TRUE) # weight the sample by inverse of cluster size

# More example with Seurat and SingelCellExperiment object as input, see vignettes.
}
