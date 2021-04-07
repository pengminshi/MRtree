#' Muti-resolution reconciled tree (MRtree)
#'
#' The method take as input the flat clusterings obtained across multiple resolution parameters. It reconcil the clusterings to produce a hierarchical cluster tree structure.
#'
#' @param x Object containing multi-resolution clustering results
#' **data source**
#' Building a reconciled hierarchical clustering tree requires information about which cluster each sample has been assinged to at different resolutions. This information can be supplied in various forms, as a matrix with rows for samples and columns for partitions, data.frame or more specialised object. In cases for data frame and objects it must contain numeric columns with the naming structure `PXS` where `P` is a prefix indicating that the column contains clustering information, `X` is a numeric value indicating the clustering resolution and `S` is any additional suffix in indicate the target columns used as input clustering labels. For `SingleCellExperiment` objects this information must be in the `colData` slot and for `Seurat` objects it must be in the `meta.data` slot.
#' @param ... Other parameters
#'
#' @return A list containing  \describe{
#'   \item{labelmat.mrtree}{The Reconciled tree saved as a label matrix, with duplicated layers omited.}
#'   \item{labelmat.recon}{The full reconciled tree save in a label matrix}
#'   \item{labelmat.flat}{The initial flat clustering (cluster tree) as input for MRtree algorithm}
#'   \item{resolutions}{The corresponding clustering resolution of the initial cluster tree}
#'   \item{paths}{The unique path in the resulted reconciled hierarchical cluster tree}
#' }
#' @export
#' @examples
#' data("clust_example")
#' # matrix as input
#' out = mrtree(clust_example$clusterings)
#' # out$labelmat.mrtree
#' # data frame with given prefix and suffix as input
#' df = as.data.frame(clust_example$clusterings); colnames(df) = paste0('K_', 1:4)
#' df$other = 1:nrow(clust_example$clusterings)  # add an additional column
#' out = mrtree(df, prefix = 'K_', suffix=NULL)
#'
#' cl = cbind(clust_example$clusterings, clust_example$clusterings)
#' # add some additional noise
#' for (i in 1:ncol(cl)){
#'     cl[sample(10),i] = sample(1:length(unique(cl[,i])), 10, replace = TRUE)
#' }
#' mrtree(cl, consensus=TRUE) # mrtree with within-resolution consensus clustering
#' mrtree(cl, augment.path=TRUE) # mrtree with augmentation of NA cluster
#' mrtree(cl, sample.weighted=TRUE) # weight the sample by inverse of cluster size
#'
#' # More example with Seurat and SingelCellExperiment object as input, see vignettes.
mrtree <- function(x, ...) {
    UseMethod("mrtree")
}


#' Muti-resolution reconciled tree (MRtree)
#'
#' \code{mrtree} with label saved in a matrix as input
#'
#' @param x matrix, label matrix where rows are samples and columns are clusterings
#' @param prefix srting indicating columns containing clustering information
#' @param suffix string at the end of column names containing clustering information
#' @param max.k the maximum resolution (number of clusters) to consider in building the tree
#' @param consensus boolean, whether to perform consensus clustering within the clusterings with the same number of clusters.
#' @param sample.weighted boolean, whether to weight the samples by the size of clusters, where higher weight is given to rare clusters. False by default.
#' @param augment.path boolean, whether to augment one NA cluster in each layer to avoid shrinking the width of the tree with clustering alternate across layers. This will dramatically increase the running time. False by default.
#' @param verbose boolean, whether to show running messages. False by default.
#' @param n.cores integer, number of cores for parallel computation. 1 core by default.
#'
#' @rdname mrtree
#' @import checkmate parallel
#' @importFrom parallel detectCores mclapply
#' @importFrom checkmate assert_true
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
mrtree.matrix <- function(x, prefix = NULL, suffix = NULL, max.k = Inf, consensus = FALSE,
    sample.weighted = FALSE, augment.path = FALSE, verbose = FALSE, n.cores = 1) {
    if (n.cores > parallel::detectCores() - 1) {
        warnings("Use ", parallel::detectCores() - 1, "cores instead!")
        n.cores = parallel::detectCores() - 1
    }
    labelmat = x

    Ks = apply(labelmat, 2, function(x) length(unique(x[!is.na(x)])))  # number of clusters per resolution

    if (is.unsorted(Ks)) {
        ord = order(Ks, decreasing = F)
        labelmat = labelmat[, ord]
        Ks = Ks[ord]
    }

    if (max.k != Inf) {
        is.effective = Ks < max.k
        labelmat = labelmat[, is.effective, drop = F]
        Ks = Ks[is.effective]
    }

    labelmat.flat = labelmat  # save the initial flat clustering results
    nk = length(Ks)  # number of initial clusterings

    if (consensus) {
        # within resolution consensus clustering
        ccs.out = consensus_clustering_within_resolution(labelmat, ks = Ks, sample.weighted = sample.weighted,
            n.cores = n.cores)
        labelmat = ccs.out$labelmat.ccs  # after consensus cluster, each layer has unique number of clusters
        labelmat.ccs = labelmat
    } else {
        labelmat.ccs = NULL
    }

    if (sample.weighted) {
        labelmat.encoded = do.call(cbind, lapply(1:ncol(labelmat), function(l) {
            label_to_membership(labelmat[, l])  # n by k
        }))

        colsum = colSums(labelmat.encoded)
        sample_weights = 1/sqrt(apply(labelmat.encoded, 1, function(x) sum(colsum[x])))  # W^{-1/2}

    } else {
        sample_weights = rep(1, nrow(labelmat))  # weights for samples
    }

    # if (any(is.na(labelmat))){ # if there are missing labels, input via nearest
    # neighbor labelmat = bnstruct::knn.impute(labelmat, k=5) labelmat.imputed =
    # labelmat } else { labelmat.imputed = NULL }

    ############### initialize
    if (verbose)
        message("initilize the tree ...")

    tree = construct_tree_from_labelmat(labelmat)  # initialize tree
    labelmat.in.paths = labelmat
    paths = unique(labelmat.in.paths)  # start from all existing paths
    bad.node = get_bad_nodeset(tree)
    if (verbose) {
        message("initial size of bad nodes:", length(bad.node))
    }

    candidate.ind = which(tree$end %in% bad.node)  # include all edges to the nodes not visited
    candidate = tree[candidate.ind, ]

    # progress bar
    message("Run MRtree:")
    total_bad = length(bad.node)
    pb <- utils::txtProgressBar(min = -total_bad - 1, max = 0, style = 3)
    lowest.layer.last = NULL


    while (length(bad.node) > 0) {
        # progress bar
        utils::setTxtProgressBar(pb, -length(bad.node))

        # top down-collection, lowest three layers
        layer = as.numeric(sapply(candidate$start, function(x) strsplit(x, split = ";")[[1]][1]))
        layer.unique.sorted = sort(unique(layer), decreasing = F)
        lowest.layer = layer.unique.sorted[1:min(2, length(layer.unique.sorted))]  # among three lowest layers
        ind.lowest.layer = which(layer %in% lowest.layer)
        candidate = candidate[ind.lowest.layer, ]

        if (augment.path) {
            paths = augment_path(paths, layers = setdiff(lowest.layer, lowest.layer.last))
            lowest.layer.last = unique(c(lowest.layer.last, lowest.layer))
        }

        # calculate the cost for candidates
        candidate$cost = unlist(parallel::mclapply(1:nrow(candidate), function(i) {
            cost(node.start = candidate$start[i], node.end = candidate$end[i], paths = paths,
                labelmat = labelmat, labelmat.in.paths = labelmat.in.paths, sample_weights = sample_weights)
        }, mc.cores = n.cores))  #, mc.cores =  parallel::detectCores()-1)

        # choose the edge with minimum cost
        ind.min = which.min(candidate$cost)
        if (length(ind.min) > 1) {
            # choose the one with max count
            ind.min = ind.min[which.max(candidate$count[ind.min])]
        }
        node.start = candidate$start[ind.min]
        node.end = candidate$en[ind.min]

        if (verbose) {
            message("Select edge-> start:", node.start, ", end:", node.end, ", cost:",
                candidate$cost[ind.min])
        }

        # remove paths that has the same end node but different start node
        if (verbose)
            message("prune path ...")
        paths = prune_paths(paths = paths, node.start = node.start, node.end = node.end)

        # only update the nodes that are affected
        if (verbose)
            message("assign sample to the path ...")
        node.end.decoded = decode(node.end)
        node.start.decoded = decode(node.start)
        node.ind.affected = which(labelmat.in.paths[, node.end.decoded$layer] ==
            node.end.decoded$label & labelmat.in.paths[, node.start.decoded$layer] !=
            node.start.decoded$label)  # only the node on the eliminated paths are affected

        if (verbose) {
            message("length(node.ind.affected)=", length(node.ind.affected))
        }

        if (length(node.ind.affected) > 0) {
            checkmate::assert_true(nrow(labelmat.in.paths) == nrow(labelmat))

            labelmat.in.paths[node.ind.affected, ] = paths[assign_samples_to_paths(labelmat = labelmat[node.ind.affected,
                , drop = F], paths = paths), , drop = F]
        }

        # update tree
        if (verbose)
            message("update the tree ...")

        tree = construct_tree_from_labelmat(labelmat.in.paths)
        bad.node = get_bad_nodeset(tree)

        if (verbose) {
            message("number of bad.node = ", length(bad.node))
        }

        candidate.ind = which((tree$end %in% bad.node))  # include all edges to the nodes not visited
        candidate = tree[candidate.ind, ]

        if (augment.path) {
            # remove the path with 0 data points assigned the path that has -1 in lowest
            # layers
            paths = rbind(unique(labelmat.in.paths), unique(paths[apply(paths[, lowest.layer,
                drop = F], 1, function(x) any(x == -1)), , drop = F]))
            paths = unique(paths)
        } else {
            paths = unique(labelmat.in.paths)
        }

        if (verbose) {
            message("number of remaining paths is ", nrow(paths))
        }
    }
    utils::setTxtProgressBar(pb, -length(bad.node))
    close(pb)

    if (consensus) {
        # map the reconciled layers with the original layers
        labelmat.recon = labelmat.in.paths[, ccs.out$k.idx]
        colnames(labelmat.recon) = colnames(labelmat.flat)
        Ks.recon = apply(labelmat.recon, 2, function(y) length(unique(y)))
    } else {
        labelmat.recon = labelmat.in.paths
        Ks.recon = apply(labelmat.recon, 2, function(y) length(unique(y)))
    }

    # unique.idx = which(!duplicated(Ks.recon[-length(Ks.recon)], MARGIN = 2))  # remove the last column in labelmat.recon
    unique.idx = which(!duplicated(Ks.recon, MARGIN = 2))  # remove the last column in labelmat.recon
    labelmat.tree = labelmat.recon[, unique.idx, drop = F]
    colnames(labelmat.tree) = paste0("K", Ks.recon[unique.idx])

    resolutions = colnames(labelmat.flat)

    return(list(labelmat.mrtree = labelmat.tree, labelmat.recon = labelmat.recon,
        labelmat.flat = labelmat.flat, labelmat.ccs = labelmat.ccs, resolutions = resolutions,
        paths = paths, params = list(prefix = prefix, suffix = suffix, consensus = consensus,
            sample.weighted = sample.weighted, augment.path=augment.path, max.k = max.k)))
}

#' Muti-resolution reconciled tree (MRtree)
#'
#' \code{mrtree} with labelmatrix save as a data frame as input
#' @param x input data.frame object
#' @param prefix srting indicating columns containing clustering information
#' @param suffix string at the end of column names containing clustering information
#' @param ... other parameters
#' @rdname mrtree
#' @importFrom checkmate assert_character
#' @export
mrtree.data.frame <- function(x, prefix=NULL, suffix = NULL, ...) {

    if (is.null(prefix) & is.null(suffix)){
        warning('Prefix and suffix both missing, using all columns')
    }

    clust_cols = rep(TRUE, ncol(x))
    if (!is.null(prefix)){
        clust_cols = clust_cols & grepl(prefix, colnames(x), fixed = TRUE)
    }

    if (!is.null(suffix)) {
        clust_cols = clust_cols & grepl(suffix, colnames(x), fixed = TRUE)
    }

    if (sum(clust_cols) < 2) {
        stop(paste("Less than two column names matched the prefix: ", prefix, "and suffix: ",
            suffix), call. = FALSE)
    }

    clusterings <- as.matrix(x[, clust_cols])
    mode(clusterings) <- "numeric"

    mrtree(clusterings, prefix = prefix, suffix = suffix, ...)
}

#' Muti-resolution reconciled tree (MRtree)
#'
#' \code{mrtree} with SingleCellExperiment object as input
#'
#' @param x SingleCellExperiment object
#' @param prefix srting indicating columns containing clustering information
#' @param suffix string at the end of column names containing clustering information
#' @param ... other parameters
#' @rdname mrtree
#' @importFrom checkmate assert_class
#' @importFrom utils installed.packages
#' @export
mrtree.SingleCellExperiment <- function(x, prefix = "sc3_", suffix = "_clusters",
    ...) {

    if (!"SingleCellExperiment" %in% utils::installed.packages()) {
        BiocManager::install("SingleCellExperiment")
        if ((!"SingleCellExperiment" %in% utils::installed.packages())) {
            stop("Please install 'SingleCellExperiment' R package.")
        }
    }

    checkmate::assert_class(x, "SingleCellExperiment")
    # checkmate::assert_character(exprs, any.missing = FALSE, len = 1)

    mrtree(data.frame(x@colData), prefix, suffix, ...)
}

#' Muti-resolution reconciled tree (MRtree)
#'
#' \code{mrtree} with Seurat object as input
#' @param x Seurat object
#' @param prefix srting indicating columns containing clustering information
#' @param suffix string at the end of column names containing clustering information
#' @param ... other parameters
#' @rdname mrtree
#' @importFrom checkmate assert_class
#' @importFrom utils installed.packages
#' @export
mrtree.Seurat <- function(x, prefix = "RNA_snn_res.", suffix = NULL, ...) {
    if (!"Seurat" %in% utils::installed.packages()) {
        install.packages("Seurat")
        if (!"Seurat" %in% utils::installed.packages()) {
            stop("Please install the Seurat R package.")
        }
    }

    checkmate::assert_class(x, "Seurat")
    # checkmate::assert_character(exprs, any.missing = FALSE)

    mrtree(x@meta.data, prefix, suffix, ...)
}


#' Construct the cluster tree from multi-resolution clustering save as a label matrix
#'
#' @param labelmat sample by m label matrix, each columns represents a clustering for a certain resolution
#' @return tree saved as an edge list containing \describe{
#'     \item{start}{label of start node}
#'     \item{end}{label of end node}
#' }
construct_tree_from_labelmat <- function(labelmat) {

    nb.layer = ncol(labelmat)

    if (nb.layer <= 1) {
        stop("Error! Can not construct tree with less than two layers")
    }

    tree = NULL
    for (l in 1:(nb.layer - 1)) {

        edgelist = table_to_edgelist(labels.start = labelmat[, l], labels.end = labelmat[,
            l + 1])

        edgelist$start = paste(l, edgelist$start, sep = ";")
        edgelist$end = paste(l + 1, edgelist$end, sep = ";")

        tree = rbind(tree, edgelist)
    }

    return(tree)
}

#' Augment the path to include on additional node as alternative to each layer
#'
#' @param paths the viable tree paths, saved in a label matrix
#' @param layers the layer index
#' @param aug.name name of the augmented cluster
#'
#' @return new paths after augmentation
augment_path <- function(paths, layers, aug.name = -1) {

    if (length(layers) == 0) {
        return(paths)
    }

    nb.layers = ncol(paths)
    nb.paths = nrow(paths)

    paths.augmented = NULL
    for (layer in layers) {
        if (layer != nb.layers & !(aug.name %in% paths[, layer])) {
            if (layer == 1) {
                aug_ind2 = apply(paths[, 2:nb.layers, drop = F], 1, function(x) any(x ==
                  -1))
                paths2 = paths[!aug_ind2, 2:nb.layers, drop = F]

                paths.augmented = rbind(path_dot_product(paths1 = matrix(aug.name,
                  nrow = 1, ncol = 1), paths2 = paths2), paths.augmented)
            } else {
                aug_ind1 = apply(paths[, 1:(layer - 1), drop = F], 1, function(x) any(x ==
                  -1))
                aug_ind2 = apply(paths[, (layer + 1):nb.layers, drop = F], 1, function(x) any(x ==
                  -1))

                paths1 = cbind(paths[!aug_ind1, 1:(layer - 1), drop = F], rep(aug.name,
                  sum(!aug_ind1)))
                paths2 = paths[!aug_ind2, (layer + 1):nb.layers, drop = F]

                paths.augmented = rbind(path_dot_product(paths1 = paths1, paths2 = paths2),
                  paths.augmented)
            }
        }
    }
    paths = unique(rbind(paths, paths.augmented))

    return(paths)
}


#' Calculate the cost generate by adding the edge = (start, end) and removing all the other conflicting edge
#'
#' @param node.start start node of the edge (out-vertex)
#' @param node.end end node of the edge (in-vertex)
#' @param paths viable path
#' @param labelmat label matrix of the samples
#' @param labelmat.in.paths assigned labels of data points to the optimimum paths
#' @param sample_weights weights per sample (clustering) used in loss function (all ones by default)
#'
#' @return a scalar representing the cost
cost <- function(node.start, node.end, paths, labelmat, labelmat.in.paths, sample_weights) {

    # prune path
    paths.new = prune_paths(paths = paths, node.start = node.start, node.end = node.end)

    # only update the nodes that are affected
    node.end.decoded = decode(node.end)
    node.ind.affected = which(labelmat.in.paths[, node.end.decoded$layer] == node.end.decoded$label)

    labelmat.in.paths.affected.old = labelmat.in.paths[node.ind.affected, , drop = F]
    labelmat.in.paths.affected.new = paths.new[assign_samples_to_paths(labelmat = labelmat[node.ind.affected,
        , drop = F], paths = paths.new), ]

    cost = sum(sample_weights[node.ind.affected] * (labelmat.in.paths.affected.old !=
        labelmat.in.paths.affected.new))

    return(cost)
}


#' Convert the start end table to an edgelist
#'
#' @param labels.start label of the start nodes
#' @param labels.end labels of the end nodes
#'
#' @return and edge list \describe{
#'     \item{start}{label of start node}
#'     \item{end}{label of end node}
#'     \item{count}{number of data points in the edge}
#' }
table_to_edgelist <- function(labels.start, labels.end) {
    tab = table(labels.start, labels.end)
    rowsum = rowSums(tab)
    colsum = colSums(tab)
    norm = matrix(rep(rowsum, ncol(tab)), ncol = ncol(tab)) + matrix(rep(colsum,
        nrow(tab)), nrow = nrow(tab), byrow = T)
    norm[norm == 0] = Inf
    tab = tab/norm
    nonzero.ind = which(tab > 0)
    edgelist = data.frame(start = rownames(tab)[row(tab)[nonzero.ind]], end = colnames(tab)[col(tab)[nonzero.ind]],
        count = tab[nonzero.ind], cost = Inf)

    return(edgelist)
}


#' parse the node name as the layer and cluster
#'
#' @param node.name name of the node
#'
#' @return A list \describe{
#'     \item{layer}{the layer in the tree the node belongs to}
#'     \item{label}{the label of the cluster in the layer the node belongs to}
#' }
decode <- function(node.name) {

    layer_label = as.numeric(unlist(strsplit(node.name, ";")))

    layer = layer_label[1]
    label = layer_label[2]

    return(list(layer = layer, label = label))
}


#' Remove the paths that include the same node.end but different node.start
#'
#' @param paths Initial paths to be pruned
#' @param node.start the label of the starting node
#' @param node.end the label of the ending node
#'
#' @return A matrix whose rows contains the remaining paths after prunning
prune_paths <- function(paths, node.start, node.end) {

    node.start.decoded = decode(node.start)
    node.end.decoded = decode(node.end)

    is.conflict = (paths[, node.start.decoded$layer] != node.start.decoded$label) &
        (paths[, node.end.decoded$layer] == node.end.decoded$label)

    if (sum(is.conflict) > 0 & node.end.decoded$layer < ncol(paths)) {
        is.selected = (paths[, node.start.decoded$layer] == node.start.decoded$label) &
            (paths[, node.end.decoded$layer] == node.end.decoded$label)
        paths.to.remove.modified = path_dot_product(paths1 = paths[is.selected, 1:node.end.decoded$layer,
            drop = F], paths2 = paths[is.conflict, (1 + node.end.decoded$layer):ncol(paths),
            drop = F])
    } else {
        paths.to.remove.modified = NULL
    }

    paths = unique(rbind(paths[!is.conflict, ], paths.to.remove.modified))

    return(paths)
}

#' outer product of path prefix and suffix
#'
#' @param paths1 first path matrix
#' @param paths2 second path matrix
#'
#' @return matrix of path outer product
path_dot_product <- function(paths1, paths2) {
    paths1 = unique(paths1)
    paths2 = unique(paths2)

    nb.paths1 = nrow(paths1)
    nb.paths2 = nrow(paths2)

    if (nb.paths1 == 0 | nb.paths2 == 0) {
        return(NULL)
    }
    path.product = cbind(paths1[rep(1:nb.paths1, nb.paths2), , drop = F], paths2[rep(1:nb.paths2,
        each = nb.paths1), , drop = F])
    return(path.product)
}


#' Assign data points to the viable path that are closest to the original clustering assignment
#'
#' @param labelmat Initial multi-resolution cluster assignment, n-by-m, corresponding to m resolutions
#' @param paths The viable paths
#' @param n.cores Number of cores to use
#'
#' @return A vector of length n containing the path index to which the data points are assigned to
#'
#' @export
assign_samples_to_paths <- function(labelmat, paths, n.cores = NULL) {

    n.layers = ncol(paths)
    paths = t(paths)

    path.labels = apply(labelmat, 1, function(label) {
        min.ind = which.min(colSums((abs(paths - label) > 0)))
        if (length(min.ind) > 1) {
            min.ind = sample(min.ind, 1)
        }
        return(min.ind)
    })

    return(path.labels)
}


#' Output the set of bad node in the given tree
#'
#' @param tree a list containg start & and pairs representing a cluster tree
#'
#' @return A vector containing the name of the bad vertex
#' @importFrom stats aggregate
get_bad_nodeset <- function(tree) {
    nb_in_edge.out = stats::aggregate(tree$start, by = list(tree$end), FUN = function(x) length(unique(x)))
    bad.node = nb_in_edge.out$Group.1[nb_in_edge.out$x > 1]

    return(bad.node)
}


#' Perform consensus clustering among clusterings of same number of clusters
#'
#' @param labelmat n-by-m sample by clustering label matrix
#' @param ks vector of length m, number of clustering for each clustering
#' @param n.cores maximum number of cores used for parallel computation
#' @param sample.weighted boolean, whether adding weights to samples
#'
#' @return a list containing \describe{
#' \item{labelmat.ccs}{consensus clustering results saved in a label matrix}
#' \item{k.idx}{idx mapping the original columns to the new columns}
#' }
#' @export
consensus_clustering_within_resolution <- function(labelmat, ks = NULL, sample.weighted = FALSE,
    n.cores = NULL) {
    if (is.null(ks)) {
        ks = apply(labelmat, 2, function(x) length(unique(x[!is.na(x)])))
    }

    unique.out = unique_ind(ks)
    k_unique = unique.out$C
    k_idx = unique.out$ic
    nk = length(k_unique)

    # if (is.null(n.cores)) { n.cores = min(nk, parallel::detectCores() - 4) }

    labelmat.ccs = do.call(cbind, lapply(1:nk, function(i) {
        k = k_unique[i]
        idx = which(k_idx == i)
        if (length(idx) > 1) {
            if (sample.weighted) {
                consensus_clustering_weighted(labelmat[, idx], k = k)
            } else {
                consensus_clustering(labelmat[, idx], k = k)
            }
        } else {
            labelmat[, idx, drop = F]
        }
    }))
    return(list(labelmat.ccs = labelmat.ccs, k.idx = k_idx))
}


#' Perform consensus clustering
#'
#' @param labelmat n-by-m sample by clustering label matrix
#' @param k number of clusters for the final consensus clustering
#'
#' @return consensus clustering result saved as a label vector (NA is labels missing for all clusterings)
#'
#' @import RSpectra
#' @importFrom stats kmeans
#' @export
consensus_clustering <- function(labelmat, k = NULL) {
    m = ncol(labelmat)
    n = nrow(labelmat)
    if (m == 1) {
        return(matrix(labelmat, ncol = 1))
    }

    if (is.null(k)) {
        # if k not provided
        ks = apply(labelmat, 2, function(x) length(unique(x[!is.na(x)])))
        k = unique(ks)
        if (length(k) > 1) {
            warning("Error, clustering applied to non-k clustering")
            k = max(ks)
        }
    }

    if (k == 1) {
        return(matrix(1, nrow = nrow(labelmat), ncol = 1))
    }

    # one-hot encode the label matrix
    labelmat.encoded = do.call(cbind, lapply(1:m, function(l) {
        label_to_membership(labelmat[, l]) # n by K
    }))

    # if there are missing values, weight rows by 1/\sqrt{#assigned clustering}
    if (any(is.na(labelmat))) {
        w = sqrt(rowSums(labelmat.encoded))
        allmissing = w == 0

        labels = rep(NA, n)
        svd.out = RSpectra::svds(labelmat.encoded[!allmissing, ]/w[!allmissing],
            k = k, nv = 0)
        u.weigted = svd.out$u * rep(svd.out$d, each = nrow(svd.out$u))
        kmeans.out = stats::kmeans(u.weigted, centers = k, nstart = 5, iter.max = 30)
        labels[!allmissing] = kmeans.out$cluster
    } else {
        svd.out = RSpectra::svds(labelmat.encoded, k = k, nv = 0)
        u.weigted = svd.out$u * rep(svd.out$d, each = nrow(svd.out$u))
        kmeans.out = stats::kmeans(u.weigted, centers = k, nstart = 5, iter.max = 30)
        labels = kmeans.out$cluster
    }

    return(labels)
}

#' weighted consensus clustering
#'
#' @param labelmat n-by-m sample by clustering label matrix
#' @param k number of clusters for the final consensus clustering
#'
#' @return consensus clustering result saved as a label vector (NA is labels missing for all clusterings)
#' @importFrom  RSpectra svds
#' @importFrom stats kmeans
#' @export
consensus_clustering_weighted <- function(labelmat, k = NULL) {
    m = ncol(labelmat)
    n = nrow(labelmat)
    if (m == 1) {
        return(matrix(labelmat, ncol = 1))
    }

    if (is.null(k)) {
        # if k not provided
        ks = apply(labelmat, 2, function(x) length(unique(x[!is.na(x)])))
        k = unique(ks)
        if (length(k) > 1) {
            warning("Error, clustering applied to non-k clustering")
            k = max(ks)
        }
    }

    if (k == 1) {
        return(matrix(1, nrow = nrow(labelmat), ncol = 1))
    }

    # one-hot encode the label matrix
    labelmat.encoded = do.call(cbind, lapply(1:m, function(l) {
        label_to_membership(labelmat[, l])
    }))

    colsum = colSums(labelmat.encoded)
    weights_sqrt = 1/sqrt(apply(labelmat.encoded, 1, function(x) sum(colsum[x])))  # W^{-1/2}

    # if there are missing values, weight rows by 1/\sqrt{#assigned clustering}
    if (any(is.na(labelmat))) {
        allmissing = rowSums(labelmat.encoded) == 0
        labels = rep(NA, n)
        svd.out = RSpectra::svds(weights_sqrt[!allmissing] * labelmat.encoded[!allmissing,
            ], k = k, nv = 0)
        u.weigted = svd.out$u * rep(svd.out$d, each = nrow(svd.out$u))
        kmeans.out = stats::kmeans(u.weigted, centers = k, nstart = 5, iter.max = 30)
        labels[!allmissing] = kmeans.out$cluster
    } else {
        svd.out = RSpectra::svds(weights_sqrt * labelmat.encoded, k = k, nv = 0)
        u.weigted = svd.out$u * rep(svd.out$d, each = nrow(svd.out$u))
        kmeans.out = stats::kmeans(u.weigted, centers = k, nstart = 5, iter.max = 30)
        labels = kmeans.out$cluster
    }

    return(labels)
}

#' Convert label vector to membership matrix
#'
#' @param labels a vector of labels from K distint classes
#' @param label.names (optional) alternative label names, used for naming columns of membership matrix
#'
#' @return an n-by-K binary membership matrix
#' @import checkmate
#' @export
label_to_membership <- function(labels, label.names = NULL) {
    if (is.null(label.names)) {
        label.names = sort(unique(labels))
    } else {
        checkmate::assert_true(all(labels %in% label.names))
    }
    k = length(label.names)
    n = length(labels)
    memb = matrix(0, nrow=n, ncol=k)
    for (i in 1:k){
        memb[labels==label.names[i],i] = 1
    }
    checkmate::assert_true(all(rowSums(memb) == 1))
    return(memb)
}
