
# setup
set.seed(1)
#             ()
#       ()          ()
#       ()        ()  ()
#     ()  ()      ()  ()

# generat the hierarchical clustering structure
clusterings = matrix(NA, nrow=200, ncol=4)
clusterings[1:200, 1] = 1
clusterings[1:100, 2] = 1
clusterings[101:200, 2] = 2
clusterings[1:100, 3] = 1
clusterings[101:150, 3] = 2
clusterings[151:200, 3] = 3
clusterings[1:50, 4] = 1
clusterings[51:100, 4] = 2
clusterings[101:150, 4] = 3
clusterings[151:200, 4] = 4

# adding some random noise
clusterings[sample(10),2] = sample(1:2, 10, replace = TRUE)
clusterings[sample(10),3] = sample(1:3, 10, replace = TRUE)
clusterings[sample(10),4] = sample(1:4, 10, replace = TRUE)

ref.labels = rep(1:4, each=50)


test_that("mrtree with matrix input", {
    out = expect_error(mrtree(clusterings), NA)
    expect_equal(nrow(clusterings), nrow(out$labelmat.mrtree))
    expect_equal(nrow(out$labelmat.recon), nrow(out$labelmat.flat))
    expect_equal(ncol(out$labelmat.recon), ncol(out$labelmat.flat))
    expect_equal(nrow(clusterings), nrow(out$labelmat.flat))
    expect_equal(ncol(clusterings), ncol(out$labelmat.flat))
})

test_that("mrtree with data frame as input", {
    df = as.data.frame(clusterings)
    colnames(df) = paste0('K_', 1:4)
    df$other = 1:nrow(clusterings)  # add an additional column

    expect_warning(mrtree(df), 'Prefix and suffix both missing')
    out = expect_error(mrtree(df, prefix = 'K_'), NA)
    expect_equal(nrow(df), nrow(out$labelmat.mrtree))
    expect_equal(nrow(out$labelmat.recon), nrow(out$labelmat.flat))
    expect_equal(ncol(out$labelmat.recon), ncol(out$labelmat.flat))
})

test_that("mrtree with within resolution consensus clustering", {
    cl = cbind(clusterings, clusterings)
    # add some additional noise
    for (i in 1:ncol(cl)){
        cl[sample(10),i] = sample(1:length(unique(cl[,i])), 10, replace = TRUE)
    }
    out = expect_error(mrtree(cl, consensus=TRUE), NA)
    expect_equal(nrow(cl), nrow(out$labelmat.ccs))
})

test_that("mrtree with augmented path", {
    out = expect_error(mrtree(clusterings, augment.path = TRUE), NA)
    expect_equal(nrow(clusterings), nrow(out$labelmat.mrtree))
    expect_equal(nrow(out$labelmat.recon), nrow(out$labelmat.flat))
    expect_equal(ncol(out$labelmat.recon), ncol(out$labelmat.flat))
    expect_equal(nrow(clusterings), nrow(out$labelmat.flat))
    expect_equal(ncol(clusterings), ncol(out$labelmat.flat))
})

test_that("mrtree with sample weighting", {
    out = expect_error(mrtree(clusterings, sample.weighted = TRUE), NA)
    expect_equal(nrow(clusterings), nrow(out$labelmat.mrtree))
    expect_equal(nrow(clusterings), nrow(out$labelmat.flat))
    expect_equal(ncol(clusterings), ncol(out$labelmat.flat))
})

