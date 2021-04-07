#' A simulated single cell RNAseq data of 500 genes-by-500 cells, generated using SymSim simulation package
#' @format A list of seven attributes
#' \describe{
#'   \item{counts}{gene-by-cell count matrix}
#'   \item{metadata}{information of cells}
#'   \item{params}{simulation parameters}
#'   \item{true_counts}{true counts without noise and dropouts}
#'   \item{gene_len}{length of genes}
#'   \item{tsne_true_counts}{tSNE plot of true counts}
#'   \item{tsne_UMI_counts}{tSNE plot of observed counts}
#' }
"data_example"


#' A array of simulated multi-resolution clustering, with columns giving the labels for each clustering.
#' @format A list of
#' \describe{
#'   \item{clusterings}{A matrix of sample by clusters, labels for each clusterings}
#'   \item{ref.labels}{A vector of reference labels of the samples.}
#' }
"clust_example"
