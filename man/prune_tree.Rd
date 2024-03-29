% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_clustering_methods.R
\name{prune_tree}
\alias{prune_tree}
\title{Prune the dendrogram to have the specified number of clusters}
\usage{
prune_tree(hc, k)
}
\arguments{
\item{hc}{a dendrogram saved in hclust object}

\item{k}{interger, level to cut the dendrogram that has this number of clusters}
}
\value{
a list containing \describe{
    \item{hc}{a hclust object of pruned dendrogram}
    \item{labels}{a vector of length n, cluster labels obtained at the specified
    tree cut}
}
}
\description{
Prune the dendrogram to have the specified number of clusters
}
