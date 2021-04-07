data("data_example", package = 'mrtree')
dat = data_example

test_that("SOUP clusterings with mrtree", {
    Ks = 5:8
    clust.out =  expect_error(sc_clustering.soup(exprs=dat$counts, type='count',
                                                 Ks=Ks, estimate.k = TRUE), NA)

    out = expect_error(mrtree(clust.out$labelmat, prefix = 'soup'), NA)
    expect_equal(nrow(clust.out$labelmat), nrow(out$labelmat.mrtree))
})

