## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- warning=F, fig.width=5, fig.height=3------------------------------------
library(mrtree)
# The data simulation can be performed using SymSim package with the following wrapper function,
# The following simulation code take some time to run (skip to load data with data("data_example"))
# tree = tree1(plot=T)$phylo.tree
# 
# symsim.param.zeisel = list(gene_effects_sd=2, sigma=0.2, scale_s=0.4, alpha_mean=0.5, alpha_sd=0.025,
#                            depth_mean=2e5,depth_sd=1.5e4, nPCR1=14, prop_hge=0.01, mean_hge=3)
# ngenes = 500
# ncells = 500
# truek = 8; min_popsize = floor(ncells/truek)
# i_minpop =  which(tree$tip.label=='Microglia')
# seed = 0
# 
# simu.out = generateDataSymSim(ncells=ncells, ngenes=ngenes, tree=tree,
#                               params=symsim.param.zeisel, plot.tsne=T,
#                               i_minpop=i_minpop, min_popsize=min_popsize,
#                               nevf = 50, n_de_evf =25, sigma=0.8, seed = seed)
# simu.out$tsne_UMI_counts
# # data_example = simu.out; usethis::use_data(data_example)
# dat = simu.out

data("data_example",package = 'mrtree') # load data
dat = data_example
# tsne plot of 500 cells
data_example$tsne_UMI_counts

## ---- fig.width=7, fig.height=4-----------------------------------------------
set.seed(1)

counts = dat$counts
metadata = dat$metadata
rownames(metadata) = dat$metadata$cellid
ref.labels = dat$metadata$type

resolutions = seq(0.1, 2, 0.1)^2

# clustering using Suerat 
seurat.out = sc_clustering.seurat(counts=counts, resolutions=resolutions, metadata=metadata, npcs=10,
                                  min.cells=0, min.features=0, scale.factor=10000, return.seurat.object=T,
                                  vars.to.regress=NULL, find.variable.features=F, verbose=F)

# initial cluster tree from Seurat flat clustering
plot_clustree(labelmat=seurat.out$seurat.clusters, prefix ='RNA_snn_res.', 
              ref.labels = ref.labels, plot.ref = F)

## ---- warning=F, message=F, fig.width=7, fig.height=4-------------------------
out = mrtree(seurat.out$obj)

plot_tree_with_piechart(labelmat=out$labelmat.mrtree, ref.labels=ref.labels, draw.pie.chart=TRUE,  
                        pie.size.offset=1, show.ref.labels = T, tip.label.dist = 10, bottom.margin=30 )

## ---- fig.width=5, fig.height=3-----------------------------------------------
amri.flat = sapply(1:ncol(out$labelmat.flat), function(i) AMRI(out$labelmat.flat[,i], ref.labels)$amri)
amri.recon = sapply(1:ncol(out$labelmat.mrtree), function(i) AMRI(out$labelmat.mrtree[,i], ref.labels)$amri)

df = rbind(data.frame(k=out$num.clust.flat, amri=amri.flat, method='Seurat flat'), 
           data.frame(k=apply(out$labelmat.mrtree, 2, FUN=function(x) length(unique(x))), amri=amri.recon, method='MRtree'))
ggplot2::ggplot(data=df, aes(x=k, y=amri, color=method)) + geom_line() + theme_bw()

## ---- fig.width=5, fig.height=3-----------------------------------------------
diff = get_index_per_layer(labelmat1=out$labelmat.flat, labelmat2=out$labelmat.recon, index.name='ari')
df = aggregate(diff, by=list(k=out$num.clust.recon), FUN=mean)

ggplot2::ggplot(data=df, aes(x=k, y=x)) + geom_line() + theme_bw() + labs(title='differece between initial clusters and MRtree')

