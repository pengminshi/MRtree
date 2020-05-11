#' Plot MRtree results as a dendrogram. If reference labels are provided, a pie chart is 
#' shown at each tree node, detailing the label proprotions.
#' 
#' @param labelmat a n by m label matrix, in ouput of \code{mrtree} function, \code{labelmat.mrtree}
#' @param ref.labels a factor or characteristic vector specifying the reference labels of n data points
#' @param factor.level.names a vector specifying the order in color the reference labels are shown
#' @param draw.pie.chart boolean, whether to draw the pie chart in the tree plot
#' @param pie.size.offset scalar offset the size of the pie chart (larger the smaller the pie chart is)
#' @param show.ref.labels boolean, whether to show the labels of major type at tree nodes and tips
#' @param tip.label.dist distance of the tip labels to the tree tips
#' @param node.label.dist distance of the node labels to the nodes
#' @param bottom.margin size of the bottom margin, need to be adjusted to show the full labels 
#' 
#' @importFrom data.tree as.Node as.phylo.Node
#' @importFrom ape Ntip Nnode
#' @importFrom tibble as_tibble
#' @importFrom tidytree full_join as.treedata
#' @importFrom scales hue_pal
#' @import ggtree
#' @export
plot_tree_with_piechart <- function(labelmat, ref.labels=NULL, factor.level.names=NULL,
                                    draw.pie.chart=TRUE,  pie.size.offset=1,
                                    show.ref.labels = T, tip.label.dist = 2, node.label.dist = 2, 
                                    bottom.margin=25 ){
  
  if(is.null(ref.labels)){
    ref.labels = rep('', nrow(labelmat))
    draw.pie.chart=F
  } else {
    ref.labels = gsub('-','_', ref.labels)
  }
  
  if(is.null(factor.level.names)){
    factor.level.names = sort(unique(ref.labels))
  } else {
    factor.level.names =  gsub('-','_', factor.level.names)
  }
  
  n = nrow(labelmat)
  p = ncol(labelmat)
  labelmat = matrix(paste(matrix(rep(colnames(labelmat),each=n), nrow = n), labelmat, sep='-'), nrow = n)
  
  # convert to data.tree data type
  df = as.data.frame(unique(labelmat), stringsAsFactors = F)
  df$pathString = apply(df, 1, function(x) paste(c('all', x), collapse='/'))
  tree.datatree = data.tree::as.Node(df)
  
  # convert to phylo tree
  tree.phylo = data.tree::as.phylo.Node(tree.datatree)
  
  # generate pie charts for each node and split
  ord = data.frame(node=1:(ape::Ntip(tree.phylo)+ape::Nnode(tree.phylo)), 
                   row.names =c(tree.phylo$tip.label, tree.phylo$node.label))
  df = data.frame(labelmat = c(labelmat), 
                  ref.labels = factor(rep(ref.labels, p), levels=factor.level.names))
  df = rbind(df, data.frame(labelmat = 'all', ref.labels))
  # calculate per type percentage
  pct = aggregate(df$ref.labels, by = list(df$labelmat), FUN = function(x) {
    t = table(x)
    t/sum(t)
  })
  
  pct = data.frame(pct$x, row.names = as.character(pct$Group.1))
  pct = transform(merge(pct,ord,by="row.names",all=TRUE), row.names=Row.names, Row.names=NULL) # use transform to remove the rownames
  
  # colnames(pct) = gsub('ExM.U','ExM-U',colnames(pct))
  
  # set the node size
  nodesize = aggregate(df$labelmat, by = list(df$labelmat), FUN = function(x) length(x))
  nodesize = data.frame(nodesize = nodesize$x/max(nodesize$x), 
                        node = ord[as.character(nodesize$Group.1),],
                        row.names = ord[as.character(nodesize$Group.1),])
  nodesize$nodesize = nodesize$nodesize^(1/5)/2/pie.size.offset  # rescale to reduce the difference
  
  # set the major label of the node and tips
  major.labels = data.frame(major.labels = factor(colnames(pct[,colnames(pct)!='node']), 
                                                  levels=factor.level.names)[apply(pct[,1:(ncol(pct)-1)], 1, which.max)], 
                            node = pct$node,
                            row.names = pct$node)
  
  # only plot the splits and leaf
  tab = table(tibble::as_tibble(tree.phylo)$parent)
  issplit = setdiff(names(tab[tab>1]),ord['all',1])
  isleaf = 1:ape::Ntip(tree.phylo)
  nodesize = nodesize[c(issplit, isleaf),]
  major.labels = major.labels[c(issplit, isleaf),]
  
  tree.plot = tidytree::full_join(tidytree::as.treedata(tree.phylo), 
                                  merge(major.labels, nodesize, by="node"),
                                  by = 'node')
  cols = scales::hue_pal()(length(levels(tree.plot@extraInfo$major.labels))) # approximate the default color
  
  if (draw.pie.chart){
    pointsize = 0.01
  } else{
    pointsize = 5
  }
  gg = ggtree::ggtree(tree.plot, size = 1) + layout_dendrogram() + xlim(bottom.margin,-110)
  
  if (show.ref.labels){
    gg = gg +
    geom_tippoint(aes(color = major.labels), size = pointsize)+
    geom_nodepoint(aes(color = major.labels), size = pointsize)+
    scale_color_manual(values = cols, drop = FALSE) + 
    geom_tiplab(aes(x = x+tip.label.dist, label = major.labels),  angle = 270, color='black') + 
    guides(colour = guide_legend(override.aes = list(size=5)))+
    labs( color = "Major type") 
  }
  
  if (draw.pie.chart){
    pies = ggtree::nodepie(pct, cols = 1:(ncol(pct)-1), color= cols[order(factor.level.names)])
    pies = pies[c(issplit, isleaf)]
    piesize = nodesize$nodesize
    gg = gg + geom_inset(pies, reverse_x = T, height = piesize, width = piesize)
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
plot_clustree <- function(labelmat, prefix = NULL, ref.labels=NULL, plot.ref=T, ...){
  library(ggraph) # needed for guide_edge_colourbar to work (bug)
  if (is.null (prefix) | is.null(labelmat) ){
    colnames(labelmat) = paste0('K',ncol(labelmat))
    prefix = 'K'
  }
    
  if (class(labelmat) != 'data.frame')
    labelmat = as.data.frame(labelmat)
  
  if (plot.ref==T & is.null(ref.labels)){
    warnings('No reference labels are provided!')
    plot.ref=F
  }
  
  if(plot.ref){
    labelmat$ref.labels = as.character(ref.labels)
    clustree::clustree(labelmat, prefix = prefix, prop_filter=0, node_colour = "ref.labels", 
                       node_colour_aggr = "getmode", node_label = "ref.labels",
                       node_label_aggr = "getmode" ,...) # cluster tree       
  } else {
    # do not plot labels 
    clustree::clustree(labelmat, prefix = prefix, prop_filter=0, ...) # cluster tree      
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


# get_combined_string <- function(v){
#     t = table(v)
#     label = NULL
#     for (l in 1:length(t)){
#         if (!is.null(label)){
#             label = paste(label,'\n')
#         }
#         label = paste0(label, names(t)[l],':',t[l])
#     }
#     return(label)
# }
# 
# 


# 
# plot_tree_labels <- function(labelmat, true.labels,
#     pie.size.offset=1,tip.label.dist = 2, node.label.dist = 2, bottom.margin=25){
# 
#     library(ggtree)
#     library(ape)
#     library(ggplot2)
# 
#     n = nrow(labelmat)
#     p = ncol(labelmat)
#     labelmat = matrix(paste(matrix(rep(colnames(labelmat),each=n), nrow = n), labelmat, sep='-'), nrow = n)
# 
#     # convert to data.tree data type
#     df = as.data.frame(unique(labelmat), stringsAsFactors = F)
#     df$pathString = apply(df, 1, function(x) paste(c('all_cells', x), collapse='/'))
#     tree.datatree = data.tree::as.Node(df)
# 
#     # convert to phylo tree
#     tree.phylo = data.tree::as.phylo.Node(tree.datatree)
# 
#     # generate pie charts for each node and split
#     ord = data.frame(node=1:(ape::Ntip(tree.phylo)+ape::Nnode(tree.phylo)), row.names =c(tree.phylo$tip.label, tree.phylo$node.label))
#     df = data.frame(labelmat = c(labelmat), true.labels = as.factor(rep(true.labels, p)))
#     df = rbind(df, data.frame(labelmat = 'all_cells', true.labels))
#     # calculate per type percentage
#     pct = aggregate(df$true.labels, by = list(df$labelmat), FUN = function(x) {
#         t = table(x)
#         t/sum(t)
#     })
#     pct = data.frame(pct$x, row.names = as.character(pct$Group.1))
#     pct = transform(merge(pct,ord,by="row.names",all=TRUE), row.names=Row.names, Row.names=NULL) # use transform to remove the rownames
# 
#     # set the node size
#     nodesize = aggregate(df$labelmat, by = list(df$labelmat), FUN = function(x) length(x))
#     nodesize = data.frame(nodesize = nodesize$x/max(nodesize$x), 
#                           node = ord[as.character(nodesize$Group.1),],
#                           row.names = ord[as.character(nodesize$Group.1),])
#     nodesize$nodesize = nodesize$nodesize^(1/5)/2/pie.size.offset  # rescale to reduce the difference
# 
#     # set the major label of the node and tips
#     major.labels = data.frame(major.labels = as.factor(colnames(pct[,colnames(pct)!='node']))[apply(pct[,1:(ncol(pct)-1)], 1, which.max)], 
#                               node = pct$node,
#                               row.names = pct$node)
# 
# 
#     # only plot the splits and leaf
#     tab = table(tibble::as_tibble(tree.phylo)$parent)
#     issplit = setdiff(names(tab[tab>1]),ord['all_cells',1])
#     isleaf = 1:ape::Ntip(tree.phylo)
#     nodesize = nodesize[c(issplit, isleaf),]
#     major.labels = major.labels[c(issplit, isleaf),]
# 
# 
#     tree.plot = tidytree::full_join(tidytree::as.treedata(tree.phylo), 
#                                     merge(major.labels, nodesize, by="node"),
#                                     by = 'node')
# 
#     cols = scales::hue_pal()(length(levels(tree.plot@extraInfo$major.labels))) # approximate the default color
# 
#    
# 
#     gg = ggtree(tree.plot, size = 1) + layout_dendrogram() + xlim(bottom.margin,-110) +
#                          geom_tippoint(aes(color = major.labels), size = 5)+
#                          geom_nodepoint(aes(color = major.labels), size = 5)+
#                          scale_color_manual(values = cols, drop = FALSE) + 
#                          geom_nodelab(aes(x = x-20*node.label.dist*nodesize),  angle = 0, color='black')+
#                          geom_tiplab(aes(x = x+tip.label.dist),  angle = 270, color='black') + 
#                          guides(colour = guide_legend(override.aes = list(size=5)))+
#                          labs( color = "Major type") 
#                   
#     gg
# }
#                          
#                          
#                          
# # plot count table
# plotContTable <- function(est_label, true_label, 
#                           true_label_order = NULL, est_label_order = NULL, short.names = NULL, 
#                           xlab = "Reference",ylab=NULL){
#   
#   if(!is.null(true_label_order)){
#       assertthat::assert_that(all(sort(unique(true_label))==sort(true_label_order)))
#   }
#   if (!is.null(est_label_order)){
#       assertthat::assert_that(all(sort(unique(est_label))==sort(est_label_order)))
#   }
#   if ("factor" %in% class(true_label)) {
#     true_label = droplevels(true_label)
#   }
#   if ("factor" %in% class(est_label)) {
#     est_label = droplevels(est_label)
#   }
#   if (is.null(short.names)) {
#     short.names = levels(factor(true_label))
#   }
#   cont.table <- table(true_label, est_label)
#   if (!is.null(true_label_order)){
#       cont.table = cont.table[true_label_order,]
#   }
#   if (!is.null(est_label_order)){
#       cont.table = cont.table[,est_label_order]
#   }
#   K <- ncol(cont.table)
#   sub.clusters <- paste0("cluster ", colnames(cont.table))
#   cont.table <- apply(as.matrix(cont.table), 2, as.integer)
#   cont.table <- data.frame(cont.table)
#   cont.table$Reference = factor(short.names, levels = short.names)
#   colnames(cont.table) <- c(sub.clusters, "Reference")
#   dat3 <- reshape2::melt(cont.table, id.var = "Reference")
#   grid.labels = as.character(dat3$value)
#   grid.labels[grid.labels == "0"] = ""
#   g <- ggplot2::ggplot(dat3, aes(Reference, variable)) + geom_tile(aes(fill = value)) + 
#     geom_text(aes(label = grid.labels), size = 4.5) + scale_fill_gradient(low = "white", 
#                                                                           high = "purple") + 
#     labs(y = ylab, x = xlab) + theme(panel.background = element_blank(), 
#                                      axis.line = element_blank(), axis.text.x = element_text(size = 13,angle = 90), 
#                                      axis.text.y = element_text(size = 13), axis.ticks = element_blank(),
#                                      axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
#                                      legend.position = "none")
#   return(g)
# }
# 
# 
# # visualize the similarity matrix
# plot_similarity <- function(A){
#   library(ggplot2)
#   A.melt = reshape2::melt(A)
#   ggplot(data = A.melt, aes(x=Var1, y=Var2, fill=value))+
#     geom_tile() + xlab('') + ylab('') #+ labs(title ='similarity matrix')
# }
# 
# 
