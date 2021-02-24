#' Plot MRtree results as a dendrogram. If reference labels are provided, a pie chart is
#' shown at each tree node, detailing the label proprotions.
#'
#' @param labelmat a n by m label matrix, in ouput of \code{mrtree} function, \code{labelmat.mrtree}
#' @param ref.labels a factor or characteristic vector specifying the reference labels of n data points
#' @param show.ref.labels boolean, whether to show the labels of major type at tree nodes and tips
#' @param label.order a vector specifying the order of labels for default colors
#' @param node.size scalar, the size of the pie chart / node
#' @param cols a vector of colors, one per each label
#' @param plot.piechart boolean, whether to draw the pie chart in the tree plot
#' @param tip.labels a vector of strings specifying the labels of tree leafs
#' @param tip.label.dist distance of the tip labels to the tree tips
#' @param show.branch.labels boolean, whether to show the branch labels for convenience of flipping branches
#' @param fip.branch a list of vectors each of size 2, indicating the branch labels to flip. Each time two branches are flipped.
#' @param legend.title string as legend title, default is an empty string
#' @param bottom.margin size of the bottom margin, need to be adjusted to show the full labels
#'
#' @importFrom data.tree as.Node as.phylo.Node
#' @importFrom ape Ntip Nnode
#' @importFrom tibble as_tibble
#' @importFrom tidytree full_join as.treedata
#' @importFrom ggtree ggtree nodepie layout_dendrogram
#' @export
plot_tree <- function(labelmat, ref.labels = NULL, show.ref.labels = T, label.order = NULL,
    node.size = 0.2, cols = NULL, plot.piechart = T, tip.labels = NULL, tip.label.dist = 2,
    show.branch.labels = F, flip.branch = NULL, legend.title = "", bottom.margin = 25) {

    if (is.null(ref.labels)) {
        ref.labels = rep("", nrow(labelmat))
        show.ref.labels = F
        draw.pie.chart = F
    } else {
        ref.labels = as.character(ref.labels)
        ref.labels = gsub("-", "_", ref.labels)
    }

    if (is.null(label.order)) {
        label.order = sort(unique(ref.labels))
    } else {
        label.order = gsub("-", "_", label.order)
        if (!all(label.order %in% ref.labels)) {
            warnings(sum(!label.order %in% ref.labels), "label name not if the reference labels!")
        }
    }

    if (plot.piechart) {
        pointsize = 0.01
    } else {
        pointsize = 5
    }

    n = nrow(labelmat)
    p = ncol(labelmat)

    # save in data.tree format
    labelmat = matrix(paste(matrix(rep(colnames(labelmat), each = n), nrow = n),
        labelmat, sep = "-"), nrow = n)
    df = as.data.frame(unique(labelmat), stringsAsFactors = F)
    df$pathString = apply(df, 1, function(x) paste(c("all", x), collapse = "/"))
    tree.datatree = data.tree::as.Node(df)

    # phylo tree for visualization
    tree.phylo = data.tree::as.phylo.Node(tree.datatree)

    # reference type per node
    ord = data.frame(node = 1:(ape::Ntip(tree.phylo) + ape::Nnode(tree.phylo)), row.names = c(tree.phylo$tip.label,
        tree.phylo$node.label))
    df = data.frame(labelmat = c(labelmat), ref.labels = rep(ref.labels, p))
    df = rbind(df, data.frame(labelmat = "all", ref.labels))

    # calculate per type percentage
    pct = aggregate(df$ref.labels, by = list(node = df$labelmat), FUN = function(x) {
        t = table(x)
        t/sum(t)
    })
    pct = data.frame(pct$x, row.names = pct$node, stringsAsFactors = F)
    pct = transform(merge(pct, ord, by = "row.names", all = TRUE), row.names = Row.names,
        Row.names = NULL)  # use transform to remove the rownames

    # set the node size
    nodesize = aggregate(df$labelmat, by = list(node = df$labelmat), FUN = function(x) length(x))
    nodesize = data.frame(nodesize = nodesize$x/max(nodesize$x), node = ord[as.character(nodesize$node),
        ], row.names = ord[as.character(nodesize$node), ])
    nodesize$nodesize = nodesize$nodesize^(1/8) * node.size  # rescale to reduce the difference

    # set the major label of the node and tips
    major.labels = data.frame(major.labels = colnames(pct[, colnames(pct) != "node"])[apply(pct[,
        1:(ncol(pct) - 1)], 1, which.max)], node = pct$node, row.names = pct$node)


    # only plot the splits and leaf
    tab = table(tibble::as_tibble(tree.phylo)$parent)
    issplit = setdiff(names(tab[tab > 1]), ord["all", 1])
    isleaf = 1:ape::Ntip(tree.phylo)
    nodesize = nodesize[c(issplit, isleaf), ]
    major.labels = major.labels[c(issplit, isleaf), ]
    major.labels$major.labels = factor(major.labels$major.labels, levels = label.order)

    # tree with label and size
    tree.plot = tidytree::full_join(tidytree::as.treedata(tree.phylo), merge(major.labels,
        nodesize, by = "node"), by = "node")

    # set the order of labels
    if (!is.null(cols)) {
        if (length(cols) != length(label.order)) {
            warnings("Number of color does not match the number of labels!")
        }
    } else {
        cols = gg_color_hue(length(label.order))  # approximate the default color
    }

    # plot the tree
    gg = ggtree::ggtree(tree.plot, size = 1) + ggtree::layout_dendrogram() + xlim(bottom.margin,
        -110)

    if (!is.null(flip.branch)) {
        for (i in 1:length(flip.branch)) {
            gg = ggtree::flip(tree_view = gg, node1 = which(gg$data$label == flip.branch[[i]][1]),
                node2 = which(gg$data$label == flip.branch[[i]][2]))
        }
    }

    if (show.ref.labels) {
        gg = gg + ggtree::geom_tippoint(aes(color = major.labels, size = nodesize),
            stroke = 0) + ggtree::geom_nodepoint(aes(color = major.labels, size = nodesize),
            stroke = 0) + scale_color_manual(values = cols, labels = label.order,
            drop = FALSE)
        if (!is.null(tip.labels)) {
            if (length(tip.labels) != sum(gg$data$isTip)) {
                stop("Error: leaf labels of different size with number of leaf: ",
                  ape::Ntip(tree.phylo), "!")
            }
            gg = gg + ggtree::geom_tiplab(aes(x = x + tip.label.dist, label = c(tip.labels[rank(gg$data$y[gg$data$isTip])],
                rep(NA, sum(!gg$data$isTip)))), angle = 270, color = "black")
        } else {
            gg = gg + ggtree::geom_tiplab(aes(x = x + tip.label.dist, label = major.labels),
                angle = 270, color = "black")
        }

        if (show.branch.labels) {
            gg = gg + ggtree::geom_nodelab(aes(x = x - 10, label = label), angle = 0,
                color = "black") + ggtree::geom_tiplab(aes(x = x - 10, label = label),
                angle = 0, color = "black")
        }
        gg = gg + guides(colour = guide_legend(override.aes = list(size = 5)), size = FALSE) +
            labs(color = legend.title)
    }

    if (plot.piechart) {
        pies = ggtree::nodepie(pct, cols = 1:(ncol(pct) - 1), color = cols[order(label.order)])
        pies = pies[c(issplit, isleaf)]
        piesize = nodesize$nodesize
        gg = gg + ggtree::geom_inset(pies, reverse_x = T, height = piesize, width = piesize)
    }

    gg
}

#' Plot MRtree results as a dendrogram. If reference labels are provided, a pie chart is
#' shown at each tree node, detailing the label proprotions.
#'
#' @param labelmat clustering results saved in a label matrix n-by-number of partitions
#' @param prefix string indicating columns containing clustering information
#' @param ref.labels reference labels to be shown at each tree node
#' @param plot.ref boolean wheather to color the tree node by the major type
#' according to reference labels
#'
#' @importFrom clustree clustree
#' @export
plot_clustree <- function(labelmat, prefix = NULL, ref.labels = NULL, plot.ref = T,
    ...) {
    library(ggraph)  # needed for guide_edge_colourbar to work (bug)
    if (is.null(prefix) | is.null(labelmat)) {
        colnames(labelmat) = paste0("K", ncol(labelmat))
        prefix = "K"
    }

    if (class(labelmat) != "data.frame")
        labelmat = as.data.frame(labelmat)

    if (plot.ref == T & is.null(ref.labels)) {
        warnings("No reference labels are provided!")
        plot.ref = F
    }

    if (plot.ref) {
        labelmat$ref.labels = as.character(ref.labels)
        clustree::clustree(labelmat, prefix = prefix, prop_filter = 0, node_colour = "ref.labels",
            node_colour_aggr = "getmode", node_label = "ref.labels", node_label_aggr = "getmode",
            ...)  # cluster tree
    } else {
        # do not plot labels
        clustree::clustree(labelmat, prefix = prefix, prop_filter = 0, ...)  # cluster tree
    }
}


#' Get the mode in the vector
#'
#' @param v a vector of numeric or character
#' @return a scalar or character representing the mode of the vector
#' @export
getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}


#' UMAP plot
#'
#' @param X n-by-p expression matrix
#' @param labels vector of sample labels
#' @param pca umap parameter, dimension of pca
#' @param n_components umap parameter, dimension of low dimensional embedding space (default 2)
#' @param n_neighbors umap parameter, number of neighbors for nearest neighbor graph
#' @param min_dist umap parameter, minimum distance of point to its nearest neighbor in the
#' embedding space.
#' @param point.size numeric scalar, point size in the plot, 0.3 by default.
#' @param alpha numeric, transparency of the points in the plot, by default alpha=1 with no transparency
#' @param title string, title of the plot, NULL by default
#' @param legend.name string, legend name, by default is 'labels'
#' @param cols vector of colors, length should the same as cardinality of labels, by default NULL
#' @param emb embedding of the UMAP if provided, by default NULL
#' @param seed random seed, by default 0
#'
#' @return A list of \describe{
#'  \item{p}{umap plot}
#'  \item{emb}{umap embedding matrix}
#' }
#'
#' @import ggplot2
#' @importFrom uwot umap
#' @export
plot_umap <- function(X = NULL, labels = NULL, pca = 50, n_components = 2, n_neighbors = 30,
    min_dist = 0.1, point.size = 0.3, alpha = 1, title = NULL, legend.name = "labels",
    cols = NULL, emb = NULL, seed = 0) {
    library(ggplot2)

    if (is.null(X) & is.null(emb)) {
        stop("data not provided!")
    }

    set.seed(seed)

    if (is.null(emb)) {
        if (!is.null(pca)) {
            if (pca > ncol(X)/2) {
                pca = NULL
            }
        }
        emb = uwot::umap(X, n_neighbors = n_neighbors, n_components = n_components,
            min_dist = min_dist, pca = pca)
    }

    df = data.frame(umap1 = emb[, 1], umap2 = emb[, 2], labels = if (!is.null(labels))
        labels else rep(0, nrow(X)))
    p = ggplot(df, aes(x = umap1, y = umap2)) + geom_point(col = "black", size = point.size,
        stroke = 0, shape = 16, alpha = alpha) + labs(x = "UMAP_1", y = "UMAP_2",
        title = title) + theme_light() + theme(plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

    if (!is.null(labels)) {
        if (is.null(legend.name)) {
            legend.name = "labels"
        }

        if (is.null(cols)) {
            cols = gg_color_hue(length(unique(labels)))
        }

        p = p + geom_point(aes(colour = labels), size = point.size, stroke = 0, shape = 16,
            alpha = alpha) + scale_color_manual(values = cols) + guides(col = guide_legend(ncol = 1,
            title = legend.name, override.aes = list(size = 5)))
    }
    list(p = p, emb = emb)
}



#' Generate colors that approximate the default ggplot colors
#'
#' @param  n integer, number of colors
#'
#' @return a vector of length n, of strings giving n colors
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


#' Plot the confusion matrix via heatmap, with true labels in columns
#'
#' @param est_label vector of estimated labels (rows)
#' @param true_label vector of actual labels (columns)
#' @param true_label_order order of column labels if provided
#' @param est_label_order order of row labels if provided
#' @param short.names short names of the true labels to be shown
#' @param xlab x axis label
#' @param ylab y axis label
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @return ggplot object
plotContTable <- function(est_label, true_label, true_label_order = NULL, est_label_order = NULL,
    short.names = NULL, xlab = "Reference", ylab = NULL) {

    library(ggplot2)
    if (!is.null(true_label_order)) {
        assertthat::assert_that(all(sort(unique(true_label)) == sort(true_label_order)))
        true_label = factor(true_label, levels = true_label_order)
    }
    if (!is.null(est_label_order)) {
        assertthat::assert_that(all(sort(unique(est_label)) == sort(est_label_order)))
        # est_label = factor(est_label, levels=est_label_order)
    }
    if (is.null(short.names)) {
        short.names = levels(factor(true_label))
    }
    cont.table <- table(true_label, est_label)
    if (!is.null(true_label_order)) {
        cont.table = cont.table[true_label_order, ]
    }
    if (!is.null(est_label_order)) {
        cont.table = cont.table[, est_label_order]
    }
    K <- ncol(cont.table)
    sub.clusters <- paste0("cluster ", colnames(cont.table))
    cont.table <- apply(as.matrix(cont.table), 2, as.integer)
    cont.table <- data.frame(cont.table)
    cont.table$Reference = factor(short.names, levels = short.names)
    colnames(cont.table) <- c(sub.clusters, "Reference")
    dat3 <- reshape2::melt(cont.table, id.var = "Reference")
    grid.labels = as.character(dat3$value)
    grid.labels[grid.labels == "0"] = ""
    g <- ggplot(dat3, aes(Reference, variable)) + geom_tile(aes(fill = value)) +
        geom_text(aes(label = grid.labels), size = 4.5) + scale_fill_gradient(low = "white",
        high = "purple") + labs(y = ylab, x = xlab) + theme(panel.background = element_blank(),
        axis.line = element_blank(), axis.text.x = element_text(size = 13, angle = 90),
        axis.text.y = element_text(size = 13), axis.ticks = element_blank(), axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18), legend.position = "none")
    return(g)
}
