
####################################################################
# simulate sc data from a hierarchical tree structure
# SymSim package
####################################################################
# generate a partition tree with 8 cell types

#' Generate tree structure 
#' 
#' @param boolean, whether to generate plot
#' 
#' @return A list containing \describe{
#'     \item{phylo.tree}{Generated tree saved in a phylo object}
#'     \item{fig}{figure that shows the tree structure}
#' }
#' 
#' @importFrom data.tree as.phylo.Node as.Node
#' @importFrom yaml yaml.load
#' @export
tree1 <- function(plot=F){

	yaml="
name: brain cell
Glia cell:
    Astrocytes:
        params: NULL
    Oligodendrocytes:
        Oligo1:
            params: NULL
        Oligo2:
            params: NULL
    Microglia:
        params: NULL
Neurons:
    Neuron1:
        Neuron1-1:
            params: NULL
        Neuron1-2:
            params: NULL
    Neuron2:
        params: NULL
    Neuron3:
        params: NULL
"
	os.list = yaml::yaml.load(yaml)
	tree = data.tree::as.Node(os.list)

	if(plot) {
		message('Plot the tree structure ..')
		p = plot(tree)
	} else {
		p = NULL
	}

	tree.phylo = data.tree::as.phylo.Node(tree)

	return(list(phylo.tree = tree.phylo, fig = p))
}



#' Generate single cell data using SymSim model
#' 
#' @param ncells scalar, number of cells
#' @param ngenes scalar, number of genes
#' @param tree a phylo object specifying the tree struction where the cell type follow
#' @param params a list containing the simulation parameters used in SymSim
#' @param plot.tsne boolean whether to generate tsne plot
#' @param ... other parameters passed to SymSim function
#' 
#' @return A list of \describe{
#'     \item{count}{ncell-by-ngene count matrix}
#'     \item{params}{simulation parameters}
#'     \item{true_counts}{ncell-by-ngene count matrix without noise and dropouts}
#'     \item{gene_len}{}
#'     \item{tsne_true_counts}{tSNE plot of the true counts}
#'     \item{tsne_UMI_counts}{tSNE plot of the simulated observed counts}
#'     \item{metadata}{Information of the cells}
#' }
#' 
#' @import checkmate
#' @importFrom ape Ntip vcv.phylo
#' @importFrom SymSim SimulateTrueCounts SimulateTrueCounts PlotTsne
#' 
#' @export
generateDataSymSim <- function(ncells, ngenes, tree, params=NULL, seed=42, plot.tsne=F, ...){
	
    library(SymSim)
	args = list(...)
	params = update_list(params, args)
	
	# check the parameters
	param.names = names(params)
	checkmate::assert_subset(param.names, choices = c('min_popsize','i_minpop','nevf','n_de_evf',
	                                                  'sigma','geffect_mean','gene_effects_sd',
	                                                  'gene_effect_prob','bimod','param_realdata',
	                                                  'scale_s','prop_hge','mean_hge','protocol',
	                                                  'alpha_mean','alpha_sd','nPCR1',
	                                                  'depth_mean','depth_sd','nbatch'))

	checkmate::assert_class(tree, 'phylo')
	
	# min_popsize
	if('min_popsize' %in% param.names){
	    checkmate::assert_numeric(params[['min_popsize']], upper=floor(ncells/ape::Ntip(tree)), lower=1)
	    min_popsize = params[['min_popsize']]
	} else {
	    min_popsize = floor(ncells/ape::Ntip(tree))
	}
	
	# i_minpop
	if('i_minpop' %in% param.names){
	    checkmate::assert_numeric(params[['i_minpop']], upper=ape::Ntip(tree), lower=1)
	    i_minpop = params[['i_minpop']]
	} else {
	    i_minpop = 1
	}
	
	# nevf
	if('nevf' %in% param.names){
	    checkmate::assert_numeric(params[['nevf']], lower=3)
	    nevf = params[['nevf']]
	} else {
	    nevf = 10
	}
	
	# n_de_evf
	if('n_de_evf' %in% param.names){
	    checkmate::assert_numeric(params[['n_de_evf']], lower=0)
	    n_de_evf = params[['n_de_evf']]
	} else {
	    n_de_evf = 5 # ? 
	}
	
	# sigma
	if('sigma' %in% param.names){
	    checkmate::assert_numeric(params[['sigma']], lower=0)
	    sigma = params[['sigma']]
	} else {
	    sigma = 0.5 
	}
	
	# geffect_mean
	if('geffect_mean' %in% param.names){
	    checkmate::assert_numeric(params[['geffect_mean']])
	    geffect_mean = params[['geffect_mean']]
	} else {
	    geffect_mean = 0 
	}
	
	# gene_effects_sd
	if('gene_effects_sd' %in% param.names){
	    checkmate::assert_numeric(params[['gene_effects_sd']], lower = 0)
	    gene_effects_sd = params[['gene_effects_sd']]
	} else {
	    gene_effects_sd = 1
	}

	# gene_effect_prob
	if('gene_effect_prob' %in% param.names){
	    checkmate::assert_numeric(params[['gene_effect_prob']], lower=0, upper=1)
	    gene_effect_prob = params[['gene_effect_prob']]
	} else {
	    gene_effect_prob = 0.3
	}

	# bimod
	if('bimod' %in% param.names){
	    checkmate::assert_numeric(params[['bimod']], lower=0, upper=1)
	    bimod = params[['bimod']]
	} else {
	    bimod = 0
	}
	
	# param_realdata
	if('param_realdata' %in% param.names){
	    checkmate::assert_choice(params[['param_realdata']], choices = c('zeisel.imputed', 'zeisel.pop4'))
	    param_realdata = params[['param_realdata']]
	} else {
	    param_realdata = 'zeisel.imputed'
	}
	
	# scale_s
	if('scale_s' %in% param.names){
	    checkmate::assert_numeric(params[['scale_s']])
	    scale_s = params[['scale_s']]
	} else {
	    scale_s = 1
	}
	
	# prop_hge
	if('prop_hge' %in% param.names){
	    checkmate::assert_numeric(params[['prop_hge']], lower=0, upper=1)
	    prop_hge = params[['prop_hge']]
	} else {
	    prop_hge = 0.015
	}
	
	# mean_hge
	if('mean_hge' %in% param.names){
	    checkmate::assert_numeric(params[['mean_hge']])
	    mean_hge = params[['mean_hge']]
	} else {
	    mean_hge = 5
	}
	
	# protocol
	if('protocol' %in% param.names){
	    checkmate::assert_choice(params[['protocol']], choices=c("nonUMI", "UMI"))
	    protocol = params[['protocol']]
	} else {
	    protocol = "UMI"
	}
	
	# alpha_mean
	if('alpha_mean' %in% param.names){
	    checkmate::assert_numeric(params[['alpha_mean']])
	    alpha_mean = params[['alpha_mean']]
	} else {
	    alpha_mean = 0.1
	}
	
	# alpha_sd
	if('alpha_sd' %in% param.names){
	    checkmate::assert_numeric(params[['alpha_sd']], lower=0)
	    alpha_sd = params[['alpha_sd']]
	} else {
	    alpha_sd = 0.002
	}
	
	# nPCR1
	if('nPCR1' %in% param.names){
	    checkmate::assert_numeric(params[['nPCR1']])
	    nPCR1 = params[['nPCR1']]
	} else {
	    nPCR1 = 16
	}
	
	# depth_mean
	if('depth_mean' %in% param.names){
	    checkmate::assert_numeric(params[['depth_mean']])
	    depth_mean = params[['depth_mean']]
	} else {
	    depth_mean = 2e5
	}
	
	# depth_sd
	if('depth_sd' %in% param.names){
	    checkmate::assert_numeric(params[['depth_sd']], lower=0)
	    depth_sd = params[['depth_sd']]
	} else {
	    depth_sd = 1.5e4
	}

	# nbatch
	if('nbatch' %in% param.names){
	    checkmate::assert_numeric(params[['nbatch']], lower=1)
	    nbatch = params[['nbatch']]
	} else {
	    nbatch = 1
	}
	
	set.seed(seed)
	
	# simulate the true counts
	true_counts.out =SymSim::SimulateTrueCounts(ncells_total=ncells, ngenes=ngenes, 
	                                            min_popsize=min_popsize, i_minpop = i_minpop, 
	                                            evf_center = 1, evf_type = "discrete", nevf = nevf,
	                                            n_de_evf = n_de_evf, vary = "s", Sigma = sigma, #impulse = F,
	                                            phyla = tree, geffect_mean = geffect_mean, gene_effects_sd = gene_effects_sd,
	                                            gene_effect_prob = gene_effect_prob, bimod = bimod,
	                                            param_realdata = param_realdata, scale_s = scale_s, 
	                                            prop_hge = prop_hge, mean_hge = mean_hge, randseed=seed, SE = F)
    

	true_counts.out$cell_meta$type = tree$tip.label[true_counts.out$cell_meta$pop]
	
	if (plot.tsne){
	    true_counts_exprs = log2(true_counts.out$counts + 1)
	    tsne_true_counts = SymSim::PlotTsne(meta=true_counts.out$cell_meta, data=true_counts_exprs, evf_type="discrete", 
	                                        n_pc=20, label='type', saving = F, plotname="true counts")[[2]]
	} else{
	    tsne_true_counts = NULL
	}

	gene_len = sample(gene_len_pool, ngenes, replace = FALSE)

	observed_counts.out = SymSim::True2ObservedCounts(true_counts=true_counts.out$counts, 
	                                                  meta_cell=true_counts.out$cell_meta, protocol=protocol,
	                                                  alpha_mean = alpha_mean, alpha_sd = alpha_sd, 
	                                                  gene_len=gene_len, # lenslope = 0.02, nbins = 20, amp_bias_limit = c(-0.2, 0.2), rate_2PCR = 0.8, 
	                                                  nPCR1 = nPCR1, # nPCR2 = 10, 	LinearAmp = F, LinearAmp_coef = 2000, SE = NULL
	                                                  depth_mean=depth_mean, depth_sd=depth_sd, nbatch = nbatch)
	
	if (plot.tsne){
	    observed_counts_exprs = log2(t(t(observed_counts.out$counts) / colSums(observed_counts.out$counts)) * 10^4 + 1)
	    tsne_UMI_counts = SymSim::PlotTsne(meta=observed_counts.out$cell_meta, data=observed_counts_exprs, 
	                                       evf_type="discrete", n_pc=20, label='type', saving = F, plotname="observed counts UMI")[[2]]
	} else{
	    tsne_UMI_counts=NULL
	}

	
	types = tree$tip.label[true_counts.out$cell_meta$pop]
	metadata = true_counts.out$cell_meta
	metadata$type = types
	rownames(metadata) = observed_counts.out$cell_meta$cellid
	
	counts = observed_counts.out$counts
	colnames(counts) = observed_counts.out$cell_meta$cellid
	rownames(counts) = paste('gene',1:ngenes,sep='-')
	
	return(list(counts = counts, # gene by cell
	            metadata = metadata,
	            params = list(min_popsize=min_popsize,
	            	i_minpop=i_minpop,
	            	nevf=nevf,
	            	n_de_evf=n_de_evf,
	            	sigma=sigma,
	            	geffect_mean=geffect_mean,
	            	gene_effects_sd=gene_effects_sd,
	            	gene_effect_prob=gene_effect_prob,
	            	bimod=bimod,
	            	param_realdata=param_realdata,
	            	scale_s=scale_s,
	            	prop_hge=prop_hge,
	            	mean_hge=mean_hge,
	            	protocol=protocol,
	            	alpha_mean=alpha_mean,
	            	alpha_sd=alpha_sd,
	            	nPCR1=nPCR1,
	                depth_mean=depth_mean,
	                depth_sd=depth_sd,
	                nbatch=nbatch),
	            true_counts=true_counts.out$counts,
	            gene_len = gene_len,
	            tsne_true_counts=tsne_true_counts,
	            tsne_UMI_counts=tsne_UMI_counts
	            ))
}

update_list <- function(list.old, list.new=NULL){
	if (is.null(list.new)){
		return (list.old)
	}

	len = length(list.old)

	names.list.new = names(list.new)
	names.list.old = names(list.old)

	comb.list = c(list.old[setdiff(names.list.old , names.list.new)], # only in old list
		list.new[intersect(names.list.new, names.list.old )], # exists in both list, take the new one
		list.new[setdiff(names.list.new, names.list.old )]) # only in new list

	return(comb.list)
}



####################################################################
# simulate k cluster singel cell data with splatter package
####################################################################
generateDataSplatter <- function(n, p, K = 1, group.prob = NULL, params = NULL, dropout = F, 
                                 de.facLoc = 2, de.prob = 0.1, path = NULL, 
                                 seed = 2019,  verbose = F){
    
    if(is.null(params)){
        if(verbose) {message('Generate splatter parameters...')}
        params <- splatter::newSplatParams()
    }
    
    if(!dropout){
        if (verbose) {message('Simulate ',n,' cells ', p,' genes from ',K,' cell type(s). SEED = ',seed )}
        dropout.type = 'none'
    } else {
        if (verbose) {message('Simulate ',n,' cells ', p,' genes with dropouts, from ',K,' cell type(s). SEED = ',seed )}
        dropout.type = 'experiment'
    }
    
    if (K==1){
        sim = splatter::splatSimulate(params,  batchCells = n, nGenes = p, de.facLoc = de.facLoc,
                                      dropout.type = dropout.type, method = 'single',seed = seed)
    } else {
        if(is.null(group.prob)){
            group.prob = rep(1/K,K)#as.numeric(gtools::rdirichlet(1, rep(2,K)))
            group.prob = group.prob/sum(group.prob)
        }
        if(is.null(path)){
            sim = splatter::splatSimulate(params,  batchCells = n, nGenes = p, dropout.type = dropout.type, de.facLoc = de.facLoc ,
                                          group.prob=group.prob, de.prob = de.prob, method = 'group',seed = seed)
        } else {
            
            sim = splatter::splatSimulate(params,  batchCells = n, nGenes = p, group.prob=group.prob, de.facLoc = de.facLoc,
                                          dropout.type = dropout.type, path.from = path, de.prob = de.prob, method = 'path',seed = seed)
        }
    }
    
    counts = t(SingleCellExperiment::counts(sim))
    
    
    if(dropout){
        d = t(SummarizedExperiment::assays(sim)$Dropout)
        if (verbose){
            message('Average dropout probability =', mean(c(d)))
        }
        
    } else {
        d = NULL
    }
    # counts = counts[,colMeans(counts>0) > 0.1]
    return(list(counts = counts, labels = sim@colData@listData$Group, dropouts = d))
}
