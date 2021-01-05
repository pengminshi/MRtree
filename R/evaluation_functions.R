########################################################
# AMRI Adjusted Multiresolution Rand Index
########################################################
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
#' @import checkmate
#'
#' @export
AMRI <- function(labels1, labels2){

    checkmate::assertTRUE(length(labels1) == length(labels2))

    labels1 = as.character(labels1)
    labels2 = as.character(labels2)
    n = length(labels1)

    # TO ADD, some preprocessing steps to remove the outliers
    k1 = length(unique(labels1))
    k2 = length(unique(labels2))

    if (k1<k2){
        diff.size = get_set_D_size(labels.low.resolution=labels1, labels.high.resolution=labels2)
    } else if (k1>k2){
        diff.size  =get_set_D_size(labels.low.resolution=labels2, labels.high.resolution=labels1)
    } else {
        # k1=k2
        diff.size = (get_set_D_size(labels.low.resolution=labels1, labels.high.resolution=labels2)+
            get_set_D_size(labels.low.resolution=labels2, labels.high.resolution=labels1))/2
    }
    mri = 1 -  diff.size/n/(n-1)* 2

    # adjust for chance under permuation model
    expected.mri = expected_MRI_perm(labels1=labels1, labels2=labels2)

    amri = (mri-expected.mri) / (1-expected.mri)

    return(list(mri=mri, amri=amri))
}

# get the size of set D
# D={pairs that are in different clusters in low resolution clustering but in the
#    same cluster in high resolution clustering}
get_set_D_size <- function(labels.low.resolution, labels.high.resolution){

    class.low = unique(labels.low.resolution) # need to filter out outliers?
    k.low = length(class.low)
    size = 0

    if ( k.low == 1){
        return(size)
    }

    for (k1 in 1:(k.low-1)){
        for (k2 in (k1+1):k.low){
            # cat('k1=',k1,', k2=',k2,'\n')
            ind1 = which(labels.low.resolution==class.low[k1])
            ind2 = which(labels.low.resolution==class.low[k2])
            ind = c(ind1, ind2)

            agg.out = aggregate(as.factor(labels.high.resolution[ind]), by=list(labels.low.resolution[ind]), FUN=table)
            if (length(agg.out$Group.1) > 1){
                size = size + sum(apply(matrix(c(agg.out[,-1]), nrow=2), 2, prod))
            }
        }
    }

    return(size)
}

# expected MRI under random permutation model
expected_MRI_perm <- function(labels1, labels2){

    count.per.cluster1 = table(labels1)
    count.per.cluster2 = table(labels2)
    n = length(labels1)

    k1 = length(unique(labels1))
    k2 = length(unique(labels2))

    prob.pairs.in.same.cluster1 = sum(sapply(count.per.cluster1, function(x) x*(x-1)/2)) /n/(n-1)*2
    prob.pairs.in.same.cluster2 = sum(sapply(count.per.cluster2, function(x) x*(x-1)/2)) /n/(n-1)*2


    if (k1<k2){
        diff.prob = (1-prob.pairs.in.same.cluster1)*prob.pairs.in.same.cluster2
    } else if (k1>k2){
        diff.prob = prob.pairs.in.same.cluster1 * (1-prob.pairs.in.same.cluster2)
    } else { # k1==k2
        diff.prob =  ((1-prob.pairs.in.same.cluster1)*prob.pairs.in.same.cluster2+
            prob.pairs.in.same.cluster1 * (1-prob.pairs.in.same.cluster2))/2
    }

    expected.mri = 1 - diff.prob

    return(expected.mri)
}


########################################################
# calculate per layer accuracy
########################################################

#' Calculate the similarity between each layer of initial cluster tree and find MRTree results
#'
#' @param labelmat1 n-by-m matrix, label matrix for cluster tree 1
#' @param labelmat2 n-by-m matrix, label matrix for cluster tree 2
#' @param index.name string indicates the type of index to evaluate, one of c('ari','hamming',
#' 'ami','amri')
#'
#' @return a numeric vector of length m, index per layer
#' @export
get_index_per_layer <- function(labelmat1, labelmat2, index.name='ari'){

    if (!is.matrix(labelmat1) & (is.matrix(labelmat2))){
        labelmat1 = rep_mat(labelmat1, ncol(labelmat2))
    } else if (!is.matrix(labelmat2) & (is.matrix(labelmat1))){
        labelmat2 = rep_mat(labelmat2, ncol(labelmat1))
    } else if (!is.matrix(labelmat1) & !is.matrix(labelmat2)){
        labelmat1 = matrix(labelmat1, ncol=1)
        labelmat2 = matrix(labelmat2, ncol=1)
    }

    checkmate::assertTRUE(ncol(labelmat1)==ncol(labelmat2))

    out = sapply(1:ncol(labelmat1), function(i){
        get_index(labels1=labelmat1[,i], labels2=labelmat2[,i], index.name = index.name)
    })

    return(out)
}

rep_mat <- function(x, k){
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
#' @importFrom aricode ARI
#' @export
get_index <- function(labels1, labels2, index.name=c('ari','hamming','ami','amri')){

    index.name = match.arg(index.name)

    if( index.name =='hamming'){
        sum(abels1 != labels2)
    } else if (index.name == 'ari'){
        aricode::ARI(labels1, labels2)
    } else if (index.name == 'ami'){
        aricode::AMI(labels1, labels2)
    } else if (index.name == 'amri'){
        AMRI(labels1=labels1, labels2=labels2)$amri
    }
}


########################################################
# Compare the similarity between two trees
########################################################

# calculate the similarity matrix from a phylo tree and sample labels
get_similarity_from_tree <- function(tree, labels){
	# check the labels are in the tree tip
	assertthat::assert_that(all(labels %in% tree$tip.label))

	labels.onehot = label_onehot(labels)
	class.dist.mat = ape::cophenetic.phylo(tree)[colnames(labels.onehot), colnames(labels.onehot)]
	class.sim.mat = 1 - class.dist.mat/max(class.dist.mat) # range [0,1]

	sim.mat = labels.onehot %*%  class.sim.mat %*% t( labels.onehot )

	return(sim.mat)
}


#
get_similarity_from_labelmat <- function(labelmat){

	d = e1071::hamming.distance(labelmat)
	sim.mat = 1-d/max(d)

	return(sim.mat)
}


# get the diff of two trees, which are represented with the similarity matrices
diff_between_similarity_mat <- function(sim.mat1, sim.mat2){
	n = nrow(sim.mat1)
	sum(abs(sim.mat1-sim.mat2))/n/n
}


# one hot code the label vector
label_onehot <- function(labels, K=NULL){

	types = unique(labels)
	n.types = length(types)

	if(is.null(K)){
		K = length(types)
	} else {
		assertthat::assert_that(n.types <= K)
	}

	n = length(labels)
	memb = matrix(0, nrow=n, ncol=K)

	for( k in 1:n.types){
		ind.k = which(labels==types[k])
		memb[ind.k,k] = 1
	}

	colnames(memb)[1:n.types] = types

	assertthat::assert_that(all(rowSums(memb)==1))

	return(memb)
}


# convert phylo_tree to labelmat
tree_to_labelmat <- function(tree, ...){
    UseMethod("tree_to_labelmat", tree)
}


tree_to_labelmat.phylo <- function(tree, Ks = NULL){

    if(is.null(Ks)){
        Ks = 2:ape::Ntip(tree)
    } else {
        Ks = Ks[Ks<=ape::Ntip(tree)]
    }
    labels = as.character(tree$sample.labels)

    labelmat = do.call(cbind,lapply(Ks, function(i){
        tip.labels = dendextend::cutree(tree,k=i)
        tip.labels[labels]
    }))
    Ks = apply(labelmat, 2, function(x) length(unique(x)))
    colnames(labelmat) = Ks
    return(labelmat)
}


tree_to_labelmat.hclust <- function(tree, Ks = NULL){

    if(is.null(Ks)){
        Ks = 2:length(tree$labels)
    } else {
        Ks = Ks[Ks<=length(tree$labels)]
    }
    labels = as.character(tree$sample.labels)

    labelmat = do.call(cbind,lapply(Ks, function(i){
        tip.labels = dendextend::cutree(tree,k=i)
        tip.labels[labels]
    }))
    Ks = apply(labelmat, 2, function(x) length(unique(x)))
    colnames(labelmat) = Ks
    return(labelmat)
}

