#' Adjusted Multi-resolution Rand Index between two clusterings
#'
#' Given two partitions on the save n data points, AMRI is used to evaluate the similarity between two
#' partitions given the two (different) resolutions. It is adapted from Adjusted Rand Index to account
#' for the difference in resolution
#'
#' @param labels1 a vector of labels representing partition1 results
#' @param labels2 a vector of labels representing partition2 results
#'
#' @return A list containing \describe{
#'     \item{mri}{scalar, calculated Mutlti-resolution Rand Index}
#'     \item{amri}{MRI adjusted for chance based on permutation model.}
#' }
#'
#' @importFrom checkmate assert_true
#'
#' @export
AMRI <- function(labels1, labels2) {

    checkmate::assert_true(length(labels1) == length(labels2))

    labels1 = as.character(labels1)
    labels2 = as.character(labels2)
    n = length(labels1)

    # TO ADD, some preprocessing steps to remove the outliers
    k1 = length(unique(labels1))
    k2 = length(unique(labels2))

    if (k1 < k2) {
        diff.size = get_set_D_size(labels.low.resolution = labels1, labels.high.resolution = labels2)
    } else if (k1 > k2) {
        diff.size = get_set_D_size(labels.low.resolution = labels2, labels.high.resolution = labels1)
    } else {
        # k1=k2
        diff.size = (get_set_D_size(labels.low.resolution = labels1, labels.high.resolution = labels2) + 
            get_set_D_size(labels.low.resolution = labels2, labels.high.resolution = labels1))/2
    }
    mri = 1 - diff.size/n/(n - 1) * 2

    # adjust for chance under permuation model
    expected.mri = expected_MRI_perm(labels1 = labels1, labels2 = labels2)

    amri = (mri - expected.mri)/(1 - expected.mri)

    return(list(mri = mri, amri = amri))
}

#' get the size of set D, {pairs that are in different clusters in low resolution clustering but in the
#'    same cluster in high resolution clustering}
#'
#' @param labels.low.resolution vector of labels at lower resolution
#' @param labels.high.resolution vector of labels at higher resolution
#'
#' @return a scalar of size of set D
get_set_D_size <- function(labels.low.resolution, labels.high.resolution) {

    class.low = unique(labels.low.resolution)  # need to filter out outliers?
    k.low = length(class.low)
    size = 0

    if (k.low == 1) {
        return(size)
    }

    for (k1 in 1:(k.low - 1)) {
        for (k2 in (k1 + 1):k.low) {
            # cat('k1=',k1,', k2=',k2,'\n')
            ind1 = which(labels.low.resolution == class.low[k1])
            ind2 = which(labels.low.resolution == class.low[k2])
            ind = c(ind1, ind2)

            agg.out = aggregate(as.factor(labels.high.resolution[ind]), by = list(labels.low.resolution[ind]), 
                FUN = table)
            if (length(agg.out$Group.1) > 1) {
                size = size + sum(apply(matrix(c(agg.out[, -1]), nrow = 2), 2, prod))
            }
        }
    }

    return(size)
}

#' expected Multi-resolution Rand Index (MRI) under random permutation model
#'
#' @param labels1 first label vector used for calculating the MRI
#' @param labels2 second label vector used for calculating the MRI
#'
#' @return scalar of expected MRI
expected_MRI_perm <- function(labels1, labels2) {

    count.per.cluster1 = table(labels1)
    count.per.cluster2 = table(labels2)
    n = length(labels1)

    k1 = length(unique(labels1))
    k2 = length(unique(labels2))

    prob.pairs.in.same.cluster1 = sum(sapply(count.per.cluster1, function(x) x * 
        (x - 1)/2))/n/(n - 1) * 2
    prob.pairs.in.same.cluster2 = sum(sapply(count.per.cluster2, function(x) x * 
        (x - 1)/2))/n/(n - 1) * 2

    
    if (k1 < k2) {
        diff.prob = (1 - prob.pairs.in.same.cluster1) * prob.pairs.in.same.cluster2
    } else if (k1 > k2) {
        diff.prob = prob.pairs.in.same.cluster1 * (1 - prob.pairs.in.same.cluster2)
    } else {
        # k1==k2
        diff.prob = ((1 - prob.pairs.in.same.cluster1) * prob.pairs.in.same.cluster2 + 
            prob.pairs.in.same.cluster1 * (1 - prob.pairs.in.same.cluster2))/2
    }

    expected.mri = 1 - diff.prob

    return(expected.mri)
}



#' Stability analysis that measures the similarity between the initial flat clusterins and reconciled tree
#' at each resolution. Similarity can be ARI or AMRI
#'
#' A line plot is generated showing the similarity across resolutions, measuring the clustering stability
#' for each resolution
#'
#' @param out MRtree output, list containing the labelmat.mrtree and labelmat.flat
#' @param index.name string, measurement of similarity, that can be 'ari' or 'amri'
#'
#' @return a dataframe showing the similarity per resolution
#' @import ggplot2
#' @export
stability_plot <- function(out, index.name = "ari") {
    if (!index.name %in% c("ari", "amri")) {
        stop("Index can only be 'ari' or 'amri'")
    }
    diff = get_index_per_layer(labelmat1 = out$labelmat.flat, labelmat2 = out$labelmat.recon, 
        index.name = index.name)
    df = aggregate(diff, by = list(k = apply(out$labelmat.flat, 2, FUN = function(x) length(unique(x)))), 
        FUN = mean)

    p = ggplot2::ggplot(data = df, aes(x = k, y = x)) + geom_line() + theme_bw() + 
        labs(x = "resolutions (K)", y = paste0(toupper(index.name), " between flat clusterings and MRtree"))
    colnames(df) = c("resolution", index.name)

    return(list(df = df, plot = p))
}

#' Calculate the similarity between each layer of initial cluster tree and find MRTree results
#'
#' @param labelmat1 n-by-m matrix, label matrix for cluster tree 1
#' @param labelmat2 n-by-m matrix, label matrix for cluster tree 2
#' @param index.name string indicates the type of index to evaluate, one of c('ari','hamming',
#' 'ami','amri')
#'
#' @return a numeric vector of length m, index per layer
#' @importFrom checkmate assert_true
#' @export
get_index_per_layer <- function(labelmat1, labelmat2, index.name = "ari") {

    if (!is.matrix(labelmat1) & (is.matrix(labelmat2))) {
        labelmat1 = rep_mat(labelmat1, ncol(labelmat2))
    } else if (!is.matrix(labelmat2) & (is.matrix(labelmat1))) {
        labelmat2 = rep_mat(labelmat2, ncol(labelmat1))
    } else if (!is.matrix(labelmat1) & !is.matrix(labelmat2)) {
        labelmat1 = matrix(labelmat1, ncol = 1)
        labelmat2 = matrix(labelmat2, ncol = 1)
    }

    checkmate::assert_true(ncol(labelmat1) == ncol(labelmat2))

    out = sapply(1:ncol(labelmat1), function(i) {
        get_index(labels1 = labelmat1[, i], labels2 = labelmat2[, i], index.name = index.name)
    })

    return(out)
}


#' repeat the vector as a matrix
#'
#' @param x input vector to repeat
#' @param k number of time to repeat
#'
#' @return a matrix of k columns
rep_mat <- function(x, k) {
    matrix(rep(x, k), ncol = k, byrow = F)
}


#' Calculate the index between two clusterings
#'
#' @param labels1 label vector of clustering 1
#' @param labels2 label vector of clustering 2
#' @param index.name string indicates the type of index to evaluate, one of c('ari','hamming',
#' 'ami','amri')
#'
#' @return a scalar, calculated index value
#'
#' @export
get_index <- function(labels1, labels2, index.name = c("ari", "hamming", "ami", "amri")) {

    index.name = match.arg(index.name)

    if (index.name == "hamming") {
        hamming.distance(labels1, labels2)
    } else if (index.name == "ari") {
        adjustedRandIndex(labels1, labels2)
    } else if (index.name == "amri") {
        AMRI(labels1 = labels1, labels2 = labels2)$amri
    }
}



#' Adjusted Rand Index
#'
#' A function to compute the adjusted rand index between two classifications.
#'
#' @param x a vector containing the labels of the first classification. Must be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param y a vector containing the labels of the second classification.
#' @return a scalar with the adjusted rand index.
#' @export
adjustedRandIndex <- function(x, y) {
    x <- as.vector(x)
    y <- as.vector(y)

    if (length(x) != length(y)) 
        stop("arguments must be vectors of the same length")

    tab <- table(x, y)

    if (all(dim(tab) == c(1, 1))) 
        return(1)

    a <- sum(choose(tab, 2))
    b <- sum(choose(rowSums(tab), 2)) - a
    c <- sum(choose(colSums(tab), 2)) - a
    d <- choose(sum(tab), 2) - a - b - c

    ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + a + c)/2 - (a + b) * 
        (a + c)/(a + b + c + d))

    return(ARI)
}


#' Hamming Distances of Vectors
#'
#' f both x and y are vectors, hamming.distance returns the Hamming distance (number of different elements) between this two vectors. If x is a matrix, the Hamming distances between the rows of x are computed and y is ignored.
#'
#' @param x a vector or matrix.
#' @param y an optional vector.
#'
#' @return distance matrix
hamming.distance <- function(x, y) {

    z <- NULL

    if (is.vector(x) && is.vector(y)) {
        z <- sum(x != y)
    } else {
        z <- matrix(0, nrow = nrow(x), ncol = nrow(x))
        for (k in 1:(nrow(x) - 1)) {
            for (l in (k + 1):nrow(x)) {
                z[k, l] <- hamming.distance(x[k, ], x[l, ])
                z[l, k] <- z[k, l]
            }
        }
        dimnames(z) <- list(dimnames(x)[[1]], dimnames(x)[[1]])
    }
    z
}



#' Calculate the pairwise similarity matrix between samples, from sample labels
#' and a phylo tree of labels.
#'
#' @param tree a phylo object giving the tree structure of labels
#' @param labels a label vector for n samples
#'
#' @return a n-by-n similarity matrix showing the similarity between samples
#' given the labels and the tree structure of labels
#' @importFrom ape cophenetic.phylo
#' @importFrom checkmate assert_true
get_similarity_from_tree <- function(tree, labels) {
    # check the labels are in the tree tip
    checkmate::assert_true(all(labels %in% tree$tip.label))

    labels.onehot = label_onehot(labels)
    class.dist.mat = ape::cophenetic.phylo(tree)[colnames(labels.onehot), colnames(labels.onehot)]
    class.sim.mat = 1 - class.dist.mat/max(class.dist.mat)  # range [0,1]

    sim.mat = labels.onehot %*% class.sim.mat %*% t(labels.onehot)

    return(sim.mat)
}


#' Get sample pairwise similarity matrix from the label matrix
#'
#' @param labelmat a n-by-m labelmatrix, n samples, and m clusterings one per column
#' @return a n-by-n similarity matrix
get_similarity_from_labelmat <- function(labelmat) {

    d = hamming.distance(labelmat)
    sim.mat = 1 - d/max(d)

    return(sim.mat)
}



#' Calculate the L1 distance between two similarity matrix.
#' Get the diff of two trees, which are represented with the similarity matrices
#'
#' @param sim.mat1 first similarity matrix
#' @param sim.mat2 second similarity matrix
#'
#' @return scalar L1 distance between the two similarity matrix
diff_between_similarity_mat <- function(sim.mat1, sim.mat2) {
    n = nrow(sim.mat1)
    diff = sum(abs(sim.mat1 - sim.mat2))/n/n

    return(diff)
}


#' one-hot code the label vector
#'
#' @param labels a vector of labels
#' @param K number of different categories, by default is the number of unique labels
#'
#' @return a n-by-K one-hot encoded binary matrix
#' @import checkmate
label_onehot <- function(labels, K = NULL) {

    types = unique(labels)
    n.types = length(types)

    if (is.null(K)) {
        K = length(types)
    } else {
        checkmate::assert_true(n.types <= K)
    }

    n = length(labels)
    memb = matrix(0, nrow = n, ncol = K)

    for (k in 1:n.types) {
        ind.k = which(labels == types[k])
        memb[ind.k, k] = 1
    }

    colnames(memb)[1:n.types] = types

    return(memb)
}


#' Convert phylo_tree to label matrix
#'
#' @param tree phylo object or hclust object
#' @param ... other parameters
tree_to_labelmat <- function(tree, ...) {
    UseMethod("tree_to_labelmat", tree)
}


#' Convert phylo tree to labelmatrix
#'
#' @param tree phylo tree
#' @param Ks number of clusters per layer
#'
#' @return a label matrix
#'
#' @importFrom ape Ntip
#' @importFrom dendextend cutree
tree_to_labelmat.phylo <- function(tree, Ks = NULL) {

    if (is.null(Ks)) {
        Ks = 2:ape::Ntip(tree)
    } else {
        Ks = Ks[Ks <= ape::Ntip(tree)]
    }
    labels = as.character(tree$sample.labels)

    labelmat = do.call(cbind, lapply(Ks, function(i) {
        tip.labels = dendextend::cutree(tree, k = i)
        tip.labels[labels]
    }))
    Ks = apply(labelmat, 2, function(x) length(unique(x)))
    colnames(labelmat) = Ks
    return(labelmat)
}

#' Convert hclust tree to label matrix
#'
#' @param tree hclust object representing a dendrogram
#' @param Ks number of clusters per layer
#'
#' @return a label matrix
#' @importFrom dendextend cutree
tree_to_labelmat.hclust <- function(tree, Ks = NULL) {

    if (is.null(Ks)) {
        Ks = 2:length(tree$labels)
    } else {
        Ks = Ks[Ks <= length(tree$labels)]
    }
    labels = as.character(tree$sample.labels)

    labelmat = do.call(cbind, lapply(Ks, function(i) {
        tip.labels = dendextend::cutree(tree, k = i)
        tip.labels[labels]
    }))
    Ks = apply(labelmat, 2, function(x) length(unique(x)))
    colnames(labelmat) = Ks
    return(labelmat)
}

