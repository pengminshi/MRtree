#' Perform Single Cell data clustering using Seurat
#'
#' @param counts n.genes-by-n.cells count matrix
#' @param resolutions vector of clustering resolution paramers (input for FindClusters)
#' @param metadata a data frame containing all the cell informations (equivalent to colData)
#' @param min.cells integer, include features detected in at least this number cells
#' @param min.features integer, include cells where at least this number features are detected
#' @param scale.factor scalar, sets the scale factor for cell-level normalization
#' @param vars.to.regress a vector of strings, variables to regress out
#' @param find.variable.features T: find highly variable features using 'mean.var.plot',
#' if >1: find variable features using 'vst' to get find.variable.features genes. F: not
#' genes selection is performed
#' @param npcs integer, number of principal component to calculate. The PCS are used for build
#' neighbor graph
#' @param seurat.graph.algorithm algorithm for modularity optimization (1 = original
#' Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm;
#' 4 = Leiden algorithm). Leiden requires the leidenalg python.
#' @param build.hierarchical.tree boolean, whether to build hierarchical tree using HAC
#' from Seurat clusters with \code{max(resolution)}.
#' @param return.seurat.object boolean, whether to return Seurat object. Save memory by setting it to \code{FALSE}.
#' @param verbose boolean, whether to print messages
#'
#' @return a list containing \describe{
#'     \item{seurat.clusters}{a data frame containing all clusteringS as columns with
#'     prefix 'RNA_snn_res.'}
#'     \item{obj}{Seurat object before clustering is performed}
#'     \item{hc.tree}{hclust object resulted from hierarchical agglomerative clustering
#'     using Seurat clusters from \code{max(resolutions)}}
#' }
#' @import checkmate
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures ScaleData RunPCA
#' FindNeighbors FindClusters BuildClusterTree VariableFeatures Tool
#' @importFrom utils install.packages
#' @export
sc_clustering.seurat <- function(counts, resolutions, metadata = NULL, min.cells = 0,
    min.features = 0, scale.factor = 10000, vars.to.regress = NULL, find.variable.features = T,
    npcs = 40, seurat.graph.algorithm = 1, build.hierarchical.tree = FALSE, return.seurat.object = FALSE,
    verbose = FALSE) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        utils::install.packages("Seurat")
        if (!requireNamespace("Seurat", quietly = TRUE)) {
            stop("Please install Seurat R package!")
        }
    }

    checkmate::assert_matrix(counts, mode = "numeric", any.missing = FALSE, col.names = "unique",
        row.names = "unique", min.cols = 2)
    checkmate::assert_numeric(resolutions, lower = 0, any.missing = FALSE)  # unique, sorted

    if (is.null(metadata)) {
        metadata = data.frame(cell.id = colnames(counts), row.names = colnames(counts))
    }
    checkmate::assert_data_frame(metadata, nrows = ncol(counts), row.names = "unique",
        null.ok = FALSE)
    checkmate::assert_subset(vars.to.regress, choices = colnames(metadata), empty.ok = TRUE)

    # create the seurat object
    obj = CreateSeuratObject(counts = counts, project = "seurat object", meta.data = metadata,
        min.cells = min.cells, min.features = min.features)

    # normalizing data
    obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = scale.factor,
        verbose = verbose)

    # get highly variable genes
    if (find.variable.features == 1) {
        obj = FindVariableFeatures(obj, selection.method = "mean.var.plot", mean.cutoff = c(0.0125,
            3), dispersion.cutoff = c(0.5, Inf), verbose = verbose)
        if (verbose)
            message("Selected ", length(VariableFeatures(object = obj)), " highly variable genes by mean.var.plot.")
    } else if (find.variable.features > 1) {
        obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = find.variable.features,
            verbose = verbose)  # select given number of genes
        if (verbose)
            message("Selected ", length(VariableFeatures(object = obj)), " highly variable genes by mean.var.plot.")
    } else {
        obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = nrow(obj),
            verbose = verbose)  # select all the genes
    }

    # using linear model to remove the effects of covariates center and scale by gene
    obj = ScaleData(obj, vars.to.regress = vars.to.regress, verbose = verbose)


    obj = RunPCA(obj, features = Seurat::VariableFeatures(object = obj), npcs = npcs,
        verbose = verbose)

    # run Seurat clustering
    obj.cluster = FindNeighbors(obj, dims = 1:npcs, verbose = verbose)
    obj.cluster = FindClusters(obj.cluster, resolution = resolutions, algorithm = seurat.graph.algorithm,
        verbose = verbose)
    seurat.clusters = obj.cluster[[]]

    if (build.hierarchical.tree) {
        obj.cluster = BuildClusterTree(object = obj.cluster, verbose = verbose)
        hc.tree = Tool(object = obj.cluster, slot = "BuildClusterTree")
        hc.tree$sample.labels = obj.cluster[["seurat_clusters"]][, 1]
        # data.tree$tip.label = cluster.major.labels
    } else {
        hc.tree = NULL
    }

    if (!return.seurat.object) {
        obj.cluster = NULL
    }

    return(list(seurat.clusters = seurat.clusters, obj = obj.cluster, hc.tree = hc.tree))
}

#' Perform Single Cell data clustering using SC3 clustering pipeline
#'
#' @param exprs n.genes-by-n.cells expression matrix
#' @param Ks vector of resolution, number of clusters
#' @param type string, type of the expression matrix, choices are 'counts', 'normcounts' and 'logcounts', and default by 'counts'
#' @param colData a dataframe containing cell informations
#' @param rowData a dataframe containing gene informations
#' @param estimate.k boolean, whether to estimate optimal number of clusters by sc3
#' @param scale.factor scalar sets the scale factor for cell-level normalization
#' @param build.hierarchical.tree boolean, whether to obtain HAC results using SC3
#' @param ... other parameters input to \code{sc3} function
#'
#' @return a list containing \describe{
#'     \item{sce}{an SingleCellExperiment object containing all the clustering results}
#'     \item{hc.tree}{hclust object resulted from hierarchical agglomerative clustering
#'     using obtained by SC3}
#' }
#'
#' @importFrom ape as.phylo
#' @importFrom checkmate assert_matrix assert_numeric
#' @importFrom utils installed.packages install.packages
#' @export
sc_clustering.sc3 <- function(exprs, Ks, type = c("counts", "normcounts"), colData = NULL,
    rowData = NULL, estimate.k = FALSE, scale.factor = 10000, build.hierarchical.tree = FALSE,
    ...) {

    if (!"SC3" %in% utils::installed.packages()) {
        BiocManager::install("SC3")
        if (!"SC3" %in% utils::installed.packages()) {
            stop("Please install SC3 r package!")
        }
    }

    checkmate::assert_matrix(exprs, mode = "numeric", any.missing = FALSE, min.cols = 2)
    checkmate::assert_numeric(Ks, lower = 2, null.ok = FALSE, any.missing = FALSE)
    type = match.arg(type)

    if (type == "counts") {
        sce = create_sce_from_counts(counts = exprs, col.data = colData, row.data = rowData,
            scale.factor = scale.factor)
    } else {
        # logcounts(sceset)
        sce = create_sce_from_normcounts(normcounts = exprs, col.data = colData,
            row.data = rowData)
    }

    sce = SC3::sc3(sce, ks = Ks, ...)  #  biology = FALSE, n_cores = NULL, gene_filter = FALSE, pct_dropout_min=-1

    if (estimate.k) {
        sce = SC3::sc3_estimate_k(sce)
    }

    if (build.hierarchical.tree) {
        # using the results from max K in Ks
        hc = sce@metadata$sc3$consensus[[length(Ks)]]$hc
        prune.out = prune_tree(hc = hc, k = max(Ks))
        hc.tree = ape::as.phylo(prune.out$hc)
        hc.tree$sample.labels = prune.out$labels
    } else {
        hc.tree = NULL
    }

    return(list(sce = sce, hc.tree = hc.tree))
}


#' Create SingleCellExperiment object from count data
#'
#' @param counts gene-by-cell count matrix
#' @param col.data a dataframe containing cell informations
#' @param row.data a dataframe containing gene informations
#' @param scale.factor scalar used for per cell sum of count normalization
#'
#' @return a SingleCellExperiment object
#' @importFrom utils installed.packages
#'
create_sce_from_counts <- function(counts, col.data, row.data = NULL, scale.factor = 10000) {

    if (!"SingleCellExperiment" %in% utils::installed.packages()) {
        BiocManager::install("SingleCellExperiment")
        if (!"SingleCellExperiment" %in% utils::installed.packages()) {
            stop("Please install SingleCellExperiment r package!")
        }
    }

    if (!"SummarizedExperiment" %in% utils::installed.packages()) {
        BiocManager::install("SummarizedExperiment")
        if (!"SummarizedExperiment" %in% utils::installed.packages()) {
            stop("Please install SummarizedExperiment r package!")
        }
    }

    if (is.null(row.data)) {
        sceset <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(counts)),
            colData = col.data)
    } else {
        sceset <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(counts)),
            colData = col.data, rowData = row.data)
    }
    # this function writes to logcounts slot exprs(sceset) <-
    # log2(calculateCPM(sceset, use_size_factors = FALSE) + 1)
    SingleCellExperiment::logcounts(sceset) = log2(t(t(counts)/colSums(counts)) *
        scale.factor + 1)
    # use gene names as feature symbols
    SummarizedExperiment::rowData(sceset)$feature_symbol = rownames(sceset)
    # remove features with duplicated names
    if (is.null(row.data)) {
        sceset <- sceset[!duplicated(SummarizedExperiment::rowData(sceset)$feature_symbol),
            ]
    }
    # QC isSpike(sceset, 'ERCC') <- grepl('^ERCC-', rownames(sceset)) sceset <-
    # calculateQCMetrics(sceset, feature_controls = list('ERCC' = isSpike(sceset,
    # 'ERCC')))
    return(sceset)
}


#' Create SingleCellExperiment object from normalized data
#'
#' @param normcounts gene-by-cell matrix of normalized expressions
#' @param col.data a dataframe containing cell informations
#' @param row.data a dataframe containing gene informations
#'
#' @return a SingleCellExperiment object
#' @importFrom SingleCellExperiment SingleCellExperiment rowData
#' @importFrom utils installed.packages
create_sce_from_normcounts <- function(normcounts, col.data, row.data = NULL) {

    if (!"SingleCellExperiment" %in% utils::installed.packages()) {
        BiocManager::install("SingleCellExperiment")
        if (!"SingleCellExperiment" %in% utils::installed.packages()) {
            stop("Please install SingleCellExperiment r package!")
        }
    }

    if (!"SummarizedExperiment" %in% utils::installed.packages()) {
        BiocManager::install("SummarizedExperiment")
        if (!"SummarizedExperiment" %in% utils::installed.packages()) {
            stop("Please install SummarizedExperiment r package!")
        }
    }

    if (is.null(row.data)) {
        sceset <- SingleCellExperiment::SingleCellExperiment(assays = list(normcounts = as.matrix(normcounts)),
            colData = col.data)
    } else {
        sceset <- SingleCellExperiment::SingleCellExperiment(assays = list(normcounts = as.matrix(normcounts)),
            colData = col.data, rowData = row.data)
    }
    SingleCellExperiment::logcounts(sceset) <- SingleCellExperiment::normcounts(sceset)
    # use gene names as feature symbols
    SummarizedExperiment::rowData(sceset)$feature_symbol <- rownames(sceset)
    # remove features with duplicated names
    if (is.null(rowData)) {
        sceset <- sceset[!duplicated(SummarizedExperiment::rowData(sceset)$feature_symbol),
            ]
    }
    # QC isSpike(sceset, 'ERCC') <- grepl('^ERCC-', rownames(sceset))
    return(sceset)
}



#' Prune the dendrogram to have the specified number of clusters
#'
#' @param hc a dendrogram saved in hclust object
#' @param k interger, level to cut the dendrogram that has this number of clusters
#'
#' @return a list containing \describe{
#'     \item{hc}{a hclust object of pruned dendrogram}
#'     \item{labels}{a vector of length n, cluster labels obtained at the specified
#'     tree cut}
#' }
prune_tree <- function(hc, k) {

    lab = -as.numeric(hc$labels)
    order = -hc$order
    n = length(hc$labels)

    for (i in 1:(n - k)) {
        lab[lab %in% hc$merge[i, ]] = i
        order[order %in% hc$merge[i, ]] = i
    }

    unique.lab = unique(lab)
    labels = -(1:k)
    names(labels) = unique.lab

    order.pruned = unique(-labels[as.character(order)])
    merge.pruned = matrix(NA, nrow = k - 1, ncol = 2)
    for (i in (n - k + 1):(n - 1)) {
        node1 = if (hc$merge[i, 1] > n - k)
            hc$merge[i, 1] - n + k else labels[as.character(hc$merge[i, 1])]
        node2 = if (hc$merge[i, 2] > n - k)
            hc$merge[i, 2] - n + k else labels[as.character(hc$merge[i, 2])]
        merge.pruned[i - n + k, ] = c(node1, node2)
    }
    height.pruned = hc$height[(n - k + 1):(n - 1)] - min(hc$height[(n - k):(n - 1)])  # seq(0,1, length.out = k)[-1]#
    labels.pruned = 1:k

    hc.pruned = list(merge = merge.pruned, height = height.pruned, order = order.pruned,
        labels = labels.pruned, method = hc$method, call = hc$call, dist.method = hc$dist.method)
    class(hc.pruned) = "hclust"

    sample.labels = -labels[as.character(lab)]
    return(list(hc = hc.pruned, labels = sample.labels))
}



#' Perform Single Cell data clustering using UMAP+kmeans
#'
#' @param exprs n.genes-by-n.cells expression matrix
#' @param Ks vector of resolution, number of clusters
#' @param type string, type of the expression matrix, choices are 'count' and 'log', and default by 'counts'
#' @param estimate.k boolean whether to estimate optimal number of clusters by ADPclust
#' @param scale.factor scalar sets the scale factor for cell-level normalization
#' @param subsample whether perform subsampling for each clustering
#' @param subsample.ratio ratio of subsampled size to the total sample size (applicable when subsample=TRUE)
#' @param n.components integer, UMAP output dimensionality, default 2
#' @param n.neighbors integer, number of neighbors to consider when build neighbor graph, default 30
#' @param column.prefix string, output column prefix, default 'tsnekmeans_'
#' @param n.cores number of cores used for parallel computation
#'
#' @return a list containing \describe{
#'     \item{labelmat}{a data frame, columns are clusterings for each resolution specified}
#'     \item{est.k}{integer, estimated number of clusters}
#' }
#'
#' @import parallel
#' @importFrom utils installed.packages
#' @export
sc_clustering.umap_kmeans <- function(exprs, Ks, type = c("count", "log"), estimate.k = FALSE,
    subsample = FALSE, subsample.ratio = 0.9, scale.factor = 10000, n.neighbors = 30,
    n.components = 2, column.prefix = "umapkmeans_", n.cores = 1) {
    if (!"uwot" %in% utils::installed.packages()) {
        install.packages("uwot")
        if (!"uwot" %in% utils::installed.packages()) {
            stop("Please install 'uwot' package")
        }
    }
    if (!"ADPclust" %in% utils::installed.packages()) {
        install.packages("ADPclust")
        if (!"ADPclust" %in% utils::installed.packages()) {
            stop("Please install 'ADPclust' package")
        }
    }

    type = match.arg(type)
    checkmate::assert_matrix(exprs, mode = "numeric", any.missing = FALSE, col.names = "unique",
        row.names = "unique", min.cols = 2)
    exprs = t(exprs)  # transposed as a ncell-by-ngenes matrix
    n = nrow(exprs)

    checkmate::assert_numeric(Ks, lower = 2, null.ok = FALSE, any.missing = FALSE)

    if (type == "count") {
        rowsum = rowSums(exprs)
        exprs = log2(exprs/rowsum * scale.factor + 1)
    }

    umap.out = uwot::umap(exprs, n_neighbors = n.neighbors, n_components = n.components,
        pca = 50)

    labels = do.call(cbind, parallel::mclapply(Ks, function(k) {
        if (subsample) {
            ind = sample(1:n, n * subsample.ratio)
        } else {
            ind = 1:n
        }

        labels.k = rep(NA, n)
        x = umap.out[ind, ]
        adp.out.k = ADPclust::adpclust(x, htype = "amise", centroids = "auto", nclust = k)
        kmeans.out = stats::kmeans(x, x[adp.out.k$centers, ], k)
        labels.k[ind] = kmeans.out$cluster

        return(labels.k)
    }, mc.cores = n.cores))
    colnames(labels) = paste0(column.prefix, "K", as.character(Ks))

    if (estimate.k) {
        adp.out = ADPclust::adpclust(umap.out, htype = "amise", centroids = "auto",
            nclust = Ks)
        est.k = adp.out$nclust
    } else {
        est.k = NULL
    }

    return(list(labelmat = as.data.frame(labels, stringsAsFactors = FALSE), est.k = est.k))
}


#' Perform Single Cell data clustering using SOUP
#'
#' @param exprs n.genes-by-n.cells expression matrix
#' @param Ks vector of resolution, number of clusters
#' @param type string, type of the expression matrix, choices are 'count' and 'log', and default by 'counts'
#' @param estimate.k boolean, whether to estimate optimal number of clusters by soup
#' @param pure.prop the proportion of pure cells, SOUP parameter
#' @param ext.prop the proportion of extreme neighbors for each cell, SOUP parameter
#' @param scale.factor scalar sets the scale factor for cell-level normalization
#' @param column.prefix string, output column prefix, default 'soup_'
#'
#' @return a list containing \describe{
#'     \item{labelmat}{a data frame, columns are clusterings for each resolution specified}
#'     \item{est.k}{integer, estimated number of clusters}
#' }
#'
#' @importFrom utils installed.packages
#' @export
sc_clustering.soup <- function(exprs, Ks, type = c("count", "log"), estimate.k = FALSE,
    pure.prop = 0.5, ext.prop = NULL, scale.factor = 10000, column.prefix = "soup_") {

    if (!"SOUP" %in% utils::installed.packages()) {
        devtools::install_github("lingxuez/SOUP")
        if (!"SOUP" %in% utils::installed.packages()) {
            stop("Please install SOUP package from Github lingxuez/SOUP")
        }
    }

    type = match.arg(type)
    checkmate::assert_matrix(exprs, mode = "numeric", any.missing = FALSE, col.names = "unique",
        row.names = "unique", min.cols = 2)
    exprs = t(exprs)  # transposed as a ncell-by-ngenes matrix
    checkmate::assert_numeric(Ks, lower = 2, null.ok = FALSE, any.missing = FALSE)

    if (type == "count") {
        exprs = log2(exprs/rowSums(exprs) * scale.factor + 1)
    }

    soup.out = SOUP::SOUP(exprs, type = "log", Ks = Ks)

    labels = do.call(cbind, soup.out$major.labels)
    colnames(labels) = paste0(column.prefix, "K", as.character(Ks))

    if (estimate.k) {
        cv.soup.out = SOUP::cvSOUP(exprs, Ks = Ks, type = "log", nfold = 10, nCV = 2,
            seeds = c(42, 42))
        est.k = cv.soup.out$K.cv
    } else {
        est.k = NULL
    }

    return(list(labelmat = as.data.frame(labels, stringsAsFactors = FALSE), est.k = est.k))
}



#' Perform Single Cell data clustering using SIMLR
#'
#' @param exprs n.genes-by-n.cells expression matrix
#' @param Ks vector of resolution, number of clusters
#' @param type string, type of the expression matrix, choices are 'count' and 'log', and default by 'counts'
#' @param estimate.k whether to estimate optimal number of clusters by SIMLR
#' @param scale.factor scalar sets the scale factor for cell-level normalization
#' @param cores.ratio ratio of total number of cores used for calculation
#' @param column.prefix string, output column prefix, default 'simlr_'
#' @return an SingleCellExperiment object containing all the clustering results
#'
#' @importFrom utils installed.packages
#' @export
sc_clustering.simlr <- function(exprs, Ks, type = c("count", "log"), estimate.k = FALSE,
    scale.factor = 10^4, cores.ratio = 0.9, column.prefix = "simlr_") {

    if (!"SIMLR" %in% utils::installed.packages()) {
        BiocManager::install("SIMLR")
        if (!"SIMLR" %in% utils::installed.packages()) {
            stop("Please install SIMLR r package")
        }
    }

    type = match.arg(type)
    checkmate::assert_matrix(exprs, mode = "numeric", any.missing = FALSE, col.names = "unique",
        row.names = "unique", min.cols = 2)

    checkmate::assert_numeric(Ks, lower = 2, null.ok = FALSE, any.missing = FALSE)

    if (type == "count") {
        colsum = colSums(exprs)
        exprs = t(log2(t(exprs)/colsum * scale.factor + 1))  # gene by cell
    }

    use.large.scale = nrow(exprs) > 1000
    labels = sapply(Ks, function(k) {
        if (use.large.scale) {
            simlr.out = SIMLR::SIMLR_Large_Scale(exprs, c = k, cores.ratio = cores.ratio,
                normalize = FALSE)
        } else {
            # small scale
            simlr.out = SIMLR::SIMLR(exprs, c = k, cores.ratio = cores.ratio, normalize = FALSE)
        }
        simlr.out$y$cluster
    })
    colnames(labels) = paste0(column.prefix, "K", as.character(Ks))

    if (estimate.k) {
        simlr.estk.out = SIMLR::SIMLR_Estimate_Number_of_Clusters(X = exprs, NUMC = Ks,
            cores.ratio = cores.ratio)
        est.k = Ks[which.min(simlr.estk.out$K1)]
    } else {
        est.k = NULL
    }

    return(list(labelmat = as.data.frame(labels, stringsAsFactors = FALSE), est.k = est.k))
}


#' Perform Single Cell data clustering using Hierarchical Agglomerative clustering (HAC)
#'
#' @param exprs n.genes-by-n.cells expression matrix
#' @param type string, type of the expression matrix, choices are 'count' and 'log', and default by 'counts
#' @param n.pcs integer, number of principal components to use for calculating similarity and performn HAC
#' @param method string, similarity method, one of 'ward.D', 'ward.D2', 'single', 'complete', 'average' (= UPGMA),
#'  'mcquitty' (= WPGMA), 'median' (= WPGMC) or 'centroid' (= UPGMC).
#' @param scale.factor scalar sets the scale factor for cell-level normalization
#' @param labels vector of string/factor, if provided, HAC is performed starting from clusters
#' @param verbose boolean, whether to print messages
#'
#' @return an hclust object containing the clustering results
#' @importFrom stats hclust prcomp dist
#' @export
sc_clustering.HAC <- function(exprs, n.pcs = 10, labels = NULL, method = "complete",
    type = c("count", "log"), scale.factor = 10000, verbose = FALSE) {

    exprs = t(exprs)

    if (type == "count") {
        exprs = log2(exprs/rowSums(exprs) * scale.factor + 1)
    }

    if (verbose) {
        message("Get first ", n.pcs, " PCs ...")
    }
    x = stats::prcomp(exprs)$x[, 1:n.pcs]

    if (!is.null(labels)) {
        labels = as.character(labels)
        message("Get ", length(unique(labels)), " cluster centers ..")
        agg.out = stats::aggregate(x, by = list(labels), FUN = mean)
        centers = agg.out[, -1]
        rownames(centers) = agg.out$Group.1

        if (verbose) {
            message("Get pairwise distance of centers ...")
        }
        d = stats::dist(centers)

        members = table(labels)[as.character(agg.out$Group.1)]
    } else {
        if (verbose) {
            message("Get pairwise distance matrix ...")
        }
        d = dist(x)
        labels = rownames(exprs)
        members = NULL
    }

    if (verbose) {
        message("HAC with ", method, " linkage ..")
    }

    hc.out = stats::hclust(d, method = method, members = members)
    hc.out$sample.labels = labels

    return(hc.out)
}
