data("data_example", package = 'mrtree')
dat = data_example
metadata = dat$metadata
rownames(metadata) = dat$metadata$cellid
ref.labels = dat$metadata$type
subset = sample(1:nrow(metadata), 100) # expedite sc3 process with fewer samples

test_that("SC3 clusterings with mrtree", {
    Ks = 5:8
    clust.out =  expect_error(sc_clustering.sc3(exprs=dat$counts[, subset],
                                                Ks=Ks, colData=metadata[subset,],
                                                gene_filter = FALSE,
                                                build.hierarchical.tree = TRUE,
                                                estimate.k = TRUE,
                                                n_cores=2), NA)
    expect_s4_class(clust.out$sce, 'SingleCellExperiment')
    expect_s3_class(clust.out$hc.tree, 'phylo')

    out = expect_error(mrtree(clust.out$sce), NA)
    expect_equal(ncol(clust.out$sce), nrow(out$labelmat.mrtree))
})

