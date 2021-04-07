
# data
data("data_example", package = 'mrtree')
dat = data_example
metadata = dat$metadata
rownames(metadata) = dat$metadata$cellid
ref.labels = dat$metadata$type

test_that("clusterings", {
    resolutions = seq(0.5, 1.4, length.out = 10)^2
    seurat.out = expect_error(sc_clustering.seurat(counts=dat$counts, resolutions=resolutions,
                                      metadata=metadata, npcs=10,
                                      return.seurat.object=TRUE,
                                      build.hierarchical.tree = TRUE,
                                      verbose=FALSE), NA)
    expect_s4_class(seurat.out$obj, 'Seurat')
    expect_s3_class(seurat.out$hc.tree, 'phylo')
})


test_that("mrtree on Seurat object / clusters", {
    resolutions = seq(0.5, 1.4, length.out = 10)^2
    seurat.out = sc_clustering.seurat(counts=dat$counts, resolutions=resolutions,
                                      metadata=metadata, npcs=10, return.seurat.object=TRUE)
    out = expect_error(mrtree(seurat.out$obj), NA)
    expect_equal(ncol(seurat.out$obj), nrow(out$labelmat.mrtree))
    out = expect_error(mrtree(seurat.out$seurat.clusters, prefix="RNA_snn_res."), NA)
    expect_equal(nrow(seurat.out$seurat.clusters), nrow(out$labelmat.mrtree))
})
