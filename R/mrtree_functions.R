#' MRtree main function
#' 
#' The multi-resolution clusters is input into \code{mrtree}, and the output is the constructed hierarchical cluster tree.
#' 
#' @param x object containing multi-resolution clustering results
#' **data source**
#' building a reconciled hierarchical clustering tree requires information about which cluster each sample has been assinged
#'  to at different resolutions. This information can be supplied in various forms, as a matrix, data.frame or more specialised
#'   object. In all cases the object provided must contain numeric columns with the naming structure `PXS` where `P` is 
#'   a prefix indicating that the column contains clustering information, `X` is a numeric value indicating the clustering
#'   resolution and `S` is any additional suffix to be removed. For `SingleCellExperiment` objects this information must 
#'   be in the `colData` slot and for `Seurat` objects it must be in the `meta.data` slot. For all objects except matrices 
#'   any additional columns can be used as aesthetics, for matrices an additional metadata data.frame can be supplied if 
#'   required.
#' @param prefix srting indicating columns containing clustering information
#' @param suffix string at the end of column names containing clustering information
#' @param max.k the maximum resolution (number of clusters) to consider in building the tree
#'
#' @return A list containing  \describe{
#'   \item{labelmat.mrtree}{The Reconciled tree saved as a label matrix, with duplicated layers omited.}
#'   \item{labelmat.recon}{The full reconciled tree save in a label matrix}
#'   \item{labelmat.flat}{The initial flat clustering (cluster tree) as input for MRtree algorithm}
#'   \item{resolutions}{The corresponding clustering resolution of the initial cluster tree}
#'   \item{paths}{The unique path in the resulted reconciled hierarchical cluster tree}
#'   \item{final.objective}{The final objective function}
#' }
#' 
#' @export
mrtree <- function(x, ...) {
    UseMethod("mrtree", x)
}



#' \code{mrtree} with label saved in a matrix as input
#' @rdname mrtree
#' @import checkmate
#' @export 
mrtree.matrix <- function(x, prefix = NULL, suffix = NULL, max.k = Inf, verbose = F) {
    
    checkmate::assert_matrix(x, any.missing = FALSE, min.cols = 2)  # col.names = 'unique', mode = 'numeric', 
    checkmate::assert_character(prefix, any.missing = FALSE, len = 1, null.ok = TRUE)
    checkmate::assert_character(suffix, any.missing = FALSE, len = 1, null.ok = TRUE)
    
    if (!is.null(colnames(x))) {
        
        if (!is.null(suffix)) {
            colnames(x) = gsub(suffix, "", colnames(x))
        }
        
        res_clean = gsub(prefix, "", colnames(x))
        resolutions = res_clean
        
    } else {
        resolutions = paste("res_", 1:ncol(x))
        colnames(x) = resolutions
    }
    
    num_clust = apply(x, 2, function(y) length(unique(y)))
    
    is.effective.k = num_clust <= max.k  # only consider the cluster with smaller than threshold number of clusters
    num_clust = num_clust[is.effective.k]
    resolutions = resolutions[is.effective.k]
    
    ord = order(num_clust, decreasing = F)
    
    labelmat = x[, is.effective.k][, ord]
    num_clust = num_clust[ord]
    weights = rep(1, length(num_clust))
    
    # initialize
    if (verbose) 
        message("initilize the tree ...")
    tree = construct_tree_from_labelmat(labelmat)  # initialize tree
    labelmat.in.paths = labelmat
    paths = unique(labelmat.in.paths)  #  # start from all existing paths
    
    bad.node = get_bad_nodeset(tree)
    if (verbose) {
        message("initial size of bad nodes:", length(bad.node))
    }
    
    candidate.ind = which((tree$end %in% bad.node))  # include all edges to the nodes not visited
    candidate = tree[candidate.ind, ]
    
    while (nrow(candidate) >= 1) {
        
        if (verbose) 
            message("\n")
        
        # add in the nodes at lower layer first
        layer = as.numeric(sapply(candidate$start, function(x) strsplit(x, split = ";")[[1]][1]))
        lowest.layer = min(layer)
        ind.lowest.layer = which(layer == lowest.layer)
        candidate = candidate[ind.lowest.layer, ]
        
        # calculate the cost for candidates
        candidate$cost = unlist(parallel::mclapply(1:nrow(candidate), function(i) {
            cost(node.start = candidate$start[i], node.end = candidate$end[i], paths = paths, 
                labelmat = labelmat, labelmat.in.paths = labelmat.in.paths, weights = weights)
        }, mc.cores = parallel::detectCores() - 1))  #, mc.cores =  parallel::detectCores()-1)
        
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
            message("length(node.ind.affected)=", length(node.ind.affected), "\n")
        }
        
        if (length(node.ind.affected) > 0) {
            checkmate::assert_true(nrow(labelmat.in.paths) == nrow(labelmat))
            
            labelmat.in.paths[node.ind.affected, ] = paths[asign_samples_to_paths(labelmat = labelmat[node.ind.affected, 
                ], paths = paths, weights.per.layer = weights), ]
        }
        
        # update tree
        if (verbose) 
            message("update the tree ...")
        tree = construct_tree_from_labelmat(labelmat.in.paths)
        # tree = update_tree(labelmat.new=labelmat.in.paths,
        # labelmat.old=labelmat.in.paths.old, node.ind.affected=node.ind.affected,
        # tree.init = tree)
        
        bad.node = get_bad_nodeset(tree)
        if (verbose) {
            message("number of bad.node = ", length(bad.node))
        }
        
        candidate.ind = which((tree$end %in% bad.node))  # include all edges to the nodes not visited
        candidate = tree[candidate.ind, ]
        
        paths = unique(labelmat.in.paths)  # remove the path with 0 data points assigned
        if (verbose) {
            message("number of remaining paths is ", nrow(paths))
        }
        
    }
    
    labelmat.recon = labelmat.in.paths
    colnames(labelmat.recon) = colnames(labelmat)
    
    num_clust_recon = apply(labelmat.recon, 2, function(y) length(unique(y)))
    unique_idx = which(!duplicated(num_clust_recon[-length(num_clust_recon)], MARGIN = 2)) # remove the last column in labelmat.recon
    labelmat.tree = labelmat.recon[, unique_idx]
    colnames(labelmat.tree) = paste0("K", num_clust_recon[unique_idx])
    
    final.objective = sum(labelmat.recon != labelmat)
    
    return(list(labelmat.mrtree = labelmat.tree, labelmat.flat = labelmat, labelmat.recon = labelmat.recon, 
        resolutions = resolutions, num.clust.flat = num_clust, num.clust.recon = num_clust_recon, 
        paths = paths, final.objective = final.objective))
}



#' \code{mrtree} with labelmatrix save as a data frame as input
#' @rdname mrtree
#' @import checkmate
#' @export 
mrtree.data.frame <- function(x, prefix, suffix = NULL, ...) {
    checkmate::assert_data_frame(x, col.names = "unique")
    checkmate::assert_character(prefix, any.missing = FALSE, len = 1)
    
    clust_cols <- grepl(prefix, colnames(x), fixed = TRUE)
    if (!is.null(suffix)) {
        clust_match_suffix <- grepl(suffix, colnames(x), fixed = TRUE)
        clust_cols = clust_cols & clust_match_suffix
    }
    
    if (sum(clust_cols) < 2) {
        stop(paste("Less than two column names matched the prefix: ", prefix, "and suffix: ", 
            suffix), call. = FALSE)
    }
    
    clusterings <- as.matrix(x[, clust_cols])
    mode(clusterings) <- "numeric"
    
    mrtree(clusterings, prefix, suffix, ...)
}




#' \code{mrtree} with SingleCellExperiment object as input
#' @rdname mrtree
#' @import checkmate
#' @export 
mrtree.SingleCellExperiment <- function(x, prefix='sc3_', suffix = "_clusters",...) {
    
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("The SingleCellExperiment package is missing, this must be",
             "installed for clustree to use SingleCellExperiment objects", 
            call. = FALSE)
    }
    
    checkmate::assert_class(x, "SingleCellExperiment")
    # checkmate::assert_character(exprs, any.missing = FALSE, len = 1)
    
    mrtree(data.frame(x@colData), prefix, suffix, ...)
}




#' \code{mrtree} with Seurat object as input
#' @rdname mrtree
#' @import checkmate
#' @export 
mrtree.Seurat <- function(x, prefix = 'RNA_snn_res.', suffix = NULL,...) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("The Seurat package is missing, this must be installed for ",
             "clustree to use seurat objects", 
            call. = FALSE)
    }
    
    warning("This interface is for the older seurat object in Seurat < 3.0.0 and ", 
        "may be deprecated in the future. You currently have Seurat v", packageVersion("Seurat"), 
        " installed. Consider installing a newer ", "version of Seurat and updating your object.", 
        call. = FALSE)
    
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



# initialize_paths <- function(Ks){ matrix(unlist(purrr::cross(lapply(Ks,
# function(k) 1:k))),ncol = length(Ks), byrow = T) }



#' Calculate the cost generate by adding the edge = (start, end) and removing all the other conflicting edge
#' 
#' @param node.start start node of the edge (out-vertex)
#' @param node.end end node of the edge (in-vertex)
#' @param paths viable path
#' @param labelmat.in.paths assigned labels of data points to the optimimum paths 
#' @param weights weights used in loss function (all ones by default)
#' 
#' @return a scalar representing the cost
cost <- function(node.start, node.end, paths, labelmat, labelmat.in.paths, weights) {
    
    # prune path
    paths.new = prune_paths(paths = paths, node.start = node.start, node.end = node.end)
    
    # only update the nodes that are affected
    node.end.decoded = decode(node.end)
    node.ind.affected = which(labelmat.in.paths[, node.end.decoded$layer] == node.end.decoded$label)
    # cat('length(node.ind.affected)=', length(node.ind.affected))
    
    labelmat.in.paths.affected.old = labelmat.in.paths[node.ind.affected, ]
    labelmat.in.paths.affected.new = paths.new[asign_samples_to_paths(labelmat = labelmat[node.ind.affected, 
        ], paths = paths.new, weights.per.layer = weights), ]
    
    cost = sum(labelmat.in.paths.affected.old != labelmat.in.paths.affected.new)
    
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
    
    # version 1 tab1 = tab/rowSums(tab) tab2 = t(t(tab)/colSums(tab)) tab =
    # matrix(unlist(Map(function(x,y) min(x, y), tab1, tab2)),nrow =nrow(tab),
    # dimnames = list(rownames(tab), colnames(tab)))
    
    # version 2
    rowsum = rowSums(tab)
    colsum = colSums(tab)
    norm = matrix(rep(rowsum, ncol(tab)), ncol = ncol(tab)) + matrix(rep(colsum, 
        nrow(tab)), nrow = nrow(tab), byrow = T)
    norm[norm == 0] = Inf
    tab = tab/norm
    
    # # version3 colsum = colSums(tab) tab = t(t(tab)/colsum)
    
    # version 4
    
    
    nonzero.ind = which(tab > 0)
    
    edgelist = data.frame(start = rownames(tab)[row(tab)[nonzero.ind]],
                          end = colnames(tab)[col(tab)[nonzero.ind]], 
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




#' Remove the paths that include the same node.end but node.start different
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
    num_conflict = sum(is.conflict)

    if (num_conflict > 0) {
        
        paths.to.remove.modified = NULL
        
        if (node.start.decoded$layer == 1) {
            # these are the path with different 'from' node but same 'to'
            paths.to.remove = matrix(paths[is.conflict, ], nrow = num_conflict)
            
            # modify these path to have the same from as the selected one
            paths.to.remove.modified = paths.to.remove
            paths.to.remove.modified[, node.start.decoded$layer] = node.start.decoded$label
            paths.to.remove.modified = unique(paths.to.remove.modified)
            
        } else {
            is.selected = (paths[, node.start.decoded$layer] == node.start.decoded$label) & 
                (paths[, node.end.decoded$layer] == node.end.decoded$label)

            path.selected = matrix(paths[is.selected, ], nrow = sum(is.selected))  # path selected, a vector
            
            # these are the path with different 'from' node but same 'to'
            paths.to.remove = matrix(paths[is.conflict, ], nrow = num_conflict) 
            have.common.parent = paths.to.remove[, node.start.decoded$layer - 1] %in% 
                path.selected[, node.start.decoded$layer - 1]
            
            if (sum(have.common.parent) > 0) {
                paths.to.remove.modified = matrix(paths.to.remove[have.common.parent, 
                  ], nrow = sum(have.common.parent))
                paths.to.remove.modified[, node.start.decoded$layer] = node.start.decoded$label
                paths.to.remove.modified = unique(paths.to.remove.modified)
                
            }
        }
    }
    
    paths = paths[!is.conflict, ]
    paths = rbind(paths, paths.to.remove.modified)
    
    return(paths)
}




#' Assign data points to the viable path that are closest to the original clustering assignment
#' 
#' @param labelmat Initial multi-resolution cluster assignment, n-by-m, corresponding to m resolutions
#' @param paths The viable paths
#' @param weights.per.layer the weights used for loss calculation (if supplied, all ones by default)
#' 
#' @return A vector of length n containing the path index to which the data points are assigned to
#' 
#' @import checkmate parallel
#' @export
asign_samples_to_paths <- function(labelmat, paths, weights.per.layer = NULL) {
    
    if (class(labelmat) != "matrix") {
        labelmat = matrix(labelmat, nrow = 1)
    }
    
    if (class(paths) != "matrix") {
        paths = matrix(paths, nrow = 1)
    }
    
    n.layers = ncol(paths)
    
    if (!is.null(weights.per.layer)) {
        checkmate::assert_true(length(weights.per.layer) == n.layers)
    } else {
        weights.per.layer = rep(1, n.layers)
    }
    
    paths = t(paths)
    
    path.labels = unlist(parallel::mclapply(1:nrow(labelmat), function(i) {
        label = labelmat[i, ]
        min.ind = which.min(colSums((abs(paths - label) > 0) * weights.per.layer))
        if (length(min.ind) > 1) {
            min.ind = sample(min.ind, 1)
        }
        min.ind
    }, mc.cores = parallel::detectCores() - 1))
    
    return(path.labels)
}




#' Output the set of bad node in the given tree
#' 
#' @param tree a list containg start & and pairs representing a cluster tree
#' 
#' @return A vector containing the name of the bad vertex
get_bad_nodeset <- function(tree) {
    nb_in_edge.out = aggregate(tree$start, by = list(tree$end), 
                               FUN = function(x) length(unique(x)))
    bad.node = nb_in_edge.out$Group.1[nb_in_edge.out$x > 1]
    
    return(bad.node)
}






