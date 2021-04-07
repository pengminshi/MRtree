data("clust_example")
clusterings = clust_example$clusterings
ref.labels = rep(c('A','B','C','D'), each=50)

out = mrtree(clusterings)

test_that("Plot the mrtree constructed tree", {
    expect_error(plot_tree(labelmat = out$labelmat.mrtree,
                           ref.labels = ref.labels,
                           plot.piechart = TRUE), NA)
    expect_error(plot_tree(labelmat = clusterings,
                           ref.labels = ref.labels,
                           plot.piechart = TRUE),
                 'Not an hierarchical tree structure')
})

test_that("Plot the mrtree constructed tree with no labels", {
    expect_error(plot_tree(labelmat=out$labelmat.mrtree,
                           ref.labels=NULL,
                           plot.piechart = FALSE), NA)
    expect_error(plot_tree(labelmat=out$labelmat.mrtree,
                           ref.labels=NULL,
                           plot.piechart = TRUE), NA)
})

test_that("Plot the mrtree constructed tree with numeric labels", {
    ref.labels = rep(1:4, each=50)
    expect_error(plot_tree(labelmat=out$labelmat.mrtree,
                           ref.labels=ref.labels,
                           plot.piechart = TRUE), NA)
})

test_that("Plot the mrtree constructed tree with mixed labels", {
    ref.labels = rep(c('A','B','C','1'), each=50)
    expect_error(plot_tree(labelmat=out$labelmat.mrtree,
                           ref.labels=ref.labels,
                           plot.piechart = TRUE), NA)
})


test_that("Plot the cluster tree", {
    expect_error(plot_clustree(labelmat=clusterings), NA)
    expect_error(plot_clustree(labelmat=clusterings,
                               ref.labels=ref.labels), NA)
})


test_that("Plot the cluster tree with repeated colnames", {
    clusterings1 = cbind(clusterings,clusterings[,4])
    colnames(clusterings1) = c('K1','K2','K3','K4',"K4")
    expect_error(plot_clustree(labelmat=clusterings1,
                               ref.labels=ref.labels), NA)
})
