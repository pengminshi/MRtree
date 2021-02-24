#' Event sample method for resolution parameters of modularity clustering.
#'
#' Sampling the resolution parameter via Event Sampling algorithm proposed in Jeub, L. G., et al. (2018). 'Multiresolution consensus clustering in networks.' Scientific reports 8(1): 1-16.
#'
#' @param A n-by-n sample adjacency matrix
#' @param n.res total number of resolution parameters to sample
#' @param gamma.min lower bound for the range of resolution parameters to sample from
#' @param gamma.max upper bound for the range of resolution parameters to sample from
#' @param subsample integer, number of sample used for deriving the resolutions via subsampling (default 200). This procedure is adopted to save computation especially when n is large.
#'
#' @return a vector of length \code{n.res}, sampled resolution parameters.
#'
#' @import checkmate
#' @export
modularity_event_sampling <- function(A, n.res, gamma.min = NULL, gamma.max = NULL,
    subsample = 200) {

    if (class(A) != "matrix") {
        A = as.matrix(A)
    }

    checkmate::assert_true(isSymmetric(A))

    if (subsample & nrow(A) > subsample) {
        subsamples = sample(1:nrow(A), subsample)
        A = A[subsamples, subsamples]
    }

    P = modularity(A)
    A = (A + t(A))/2
    P = A - P

    diag(A) = 0
    diag(P) = 0

    # get the discrete events where the interaction change sign
    gamma_et = div_0(A, P)
    unique.out = unique_ind(gamma_et)
    g_sample = unique.out$C
    ind = unique.out$ic

    # for each unique gamma changepoint, get corresponding beta
    m = length(g_sample)
    Asum = sum(A)
    Psum = sum(P)

    Ap = rep(0, m)
    Pp = rep(0, m)  # sum of P for entries in E+
    for (i in seq(m, 1)) {
        Ap[i] = sum(A[ind == i])
        Pp[i] = sum(P[ind == i])
        if (i != m) {
            Ap[i] = Ap[i] + Ap[i + 1]
            Pp[i] = Pp[i] + Pp[i + 1]
        }
    }

    b_sample = (g_sample * (Psum - Pp) - (Asum - Ap))/(g_sample * (Psum - 2 * Pp) +
        (2 * Ap - Asum))

    # remove the gamma change point that corresponds to the same beta
    unique.idx = which(!duplicated(b_sample))
    b_sample = b_sample[unique.idx]
    g_sample = g_sample[unique.idx]
    Pp = Pp[unique.idx]
    Ap = Ap[unique.idx]

    bmin = (gamma.min * grid_interpolant_next(gamma.min, g_sample, Psum - Pp) - grid_interpolant_next(gamma.min,
        g_sample, Asum - Ap))/(gamma.min * grid_interpolant_next(gamma.min, g_sample,
        Psum - 2 * Pp) + grid_interpolant_next(gamma.min, g_sample, 2 * Ap - Asum))
    bmax = (gamma.max * grid_interpolant_next(gamma.max, g_sample, Psum - Pp) - grid_interpolant_next(gamma.max,
        g_sample, Asum - Ap))/(gamma.max * grid_interpolant_next(gamma.max, g_sample,
        Psum - 2 * Pp) + grid_interpolant_next(gamma.max, g_sample, 2 * Ap - Asum))
    bgrid = seq(from = bmin, to = bmax, length.out = n.res)

    Pminus = grid_interpolant_next(grid = bgrid, x = b_sample, y = Psum - Pp)
    Pplus = grid_interpolant_next(grid = bgrid, x = b_sample, y = Pp)
    Aminus = grid_interpolant_next(grid = bgrid, x = b_sample, y = Asum - Ap)
    Aplus = grid_interpolant_next(grid = bgrid, x = b_sample, y = Ap)

    gammas = (Aminus + bgrid * (Aplus - Aminus))/((1 - bgrid) * Pminus + bgrid *
        Pplus)

    return(gammas)
}


#' Linear sampling.
#'
#' Sample resolutions parameters linearly withing a range
#'
#' @param n.res integer, number of resolutions to sample
#' @param res.min scalar, lower bound
#' @param res.max scalar, upper bound
#'
#' @return a vector of sampled parameters
#' @export
linear_sampling <- function(n.res, res.min = 0.001, res.max = 3) {
    return(seq(from = res.min, to = res.max, length.out = n.res))
}


#' Exponential sampling
#'
#' Sample resolutions in exponential scale (sparser for larger values)
#'
#' @param n.res integer, number of resolutions to sample
#' @param res.min scalar, lower bound
#' @param res.max scalar, upper bound
#' @param base base of the exponential scale (default 2)
#'
#' @return a vector of sampled parameters
#' @export
exponential_sampling <- function(n.res, res.min = 0.001, res.max = 3, base = 2) {
    if (res.min <= 0 | res.max < 0) {
        stop("Error: res.min should > 0!")
    }

    res_log = seq(from = log(res.min, base = base), to = log(res.max, base = base),
        length.out = n.res)
    return(base^res_log)
}



#' pointwise division of matrices such that 0/0=0
div_0 <- function(A, B) {

    Btmp = B
    Btmp[Btmp == 0] = 1

    return(A/Btmp)
}

#' Get the unique values in the matrix
#'
#' @param matrix A
#'
#' @return a list of \describe{
#' \item{C}{sorted unique values (increaing order)}
#' \item{ia}{a vector of index such that \code{C = c(A)[ia]}}
#' \item{ic}{index vector such that \code{c(A) = C[ic]}}
#' }
unique_ind <- function(A) {
    A = c(A)
    ia_unordered = which(!duplicated(A))
    C_unordered = A[ia_unordered]

    ord = order(C_unordered, decreasing = F)
    ia = ia_unordered[ord]
    C = C_unordered[ord]

    # get the ic
    ic = match(A, C)

    return(list(C = C, ia = ia, ic = ic))
}

#' Grid intepolation using the next value
#'
#' @param grid a vector of new xs for which the ys will be generated by intepolation
#' @param x a vector of reference x where y is given
#' @param y a vector of value corresponding to x
#'
#' @return a vector of intepolated y values for the grid of x
grid_interpolant_next <- function(grid, x, y) {
    if (is.unsorted(x)) {
        ord = order(x)
        x = x[ord]
        y = y[ord]
    }

    m = length(grid)
    n = length(x)

    interplant = rep(NA, m)

    interplant[grid <= x[1]] = y[1]
    interplant[grid >= x[n]] = y[n]

    grid.idx = sum(grid <= x[1]) + 1
    x.idx = grid.idx

    while (grid.idx <= m - sum(grid >= x[n])) {
        while (x[x.idx] < grid[grid.idx]) {
            x.idx = x.idx + 1
        }
        interplant[grid.idx] = y[x.idx]
        grid.idx = grid.idx + 1
    }
    return(interplant)
}

#' Modularity function given the adjacency matrix
#'
#' @param A adjacency matrix
#' @param gamma scalar, parameter for modularity calculation
#'
#' @return numeric, calculated modularity value
modularity <- function(A, gamma = 1) {
    k = matrix(apply(A, 2, sum), ncol = 1)
    d = matrix(apply(A, 1, sum), nrow = 1)
    twom = sum(k)

    B = (A + t(A))/2 - gamma/2 * (k %*% d + t(d) %*% t(k))/twom
    return(B)
}



#' get the nearest neighbor graph using Seurat
#'
#' @param counts ncell by ngene count matrix
#' @param metadata a dataframe of meta data of cells
#' @param npc principal components for constructing the nearest neighbor graph
#'
#' @return nearest neighbor graph saved as an adjacency matrix
#' @export
seurat_get_nn_graph <- function(counts, metadata = NULL, npcs = 10, ...) {

    require(Seurat, quietly = T)

    min.cells = 0
    min.features = 0
    scale.factor = 10^4
    find.variable.features = F
    vars.to.regress = NULL
    verbose = F

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


    obj = RunPCA(obj, features = VariableFeatures(object = obj), npcs = npcs, verbose = verbose)

    # run Seurat clustering
    obj.cluster = FindNeighbors(obj, dims = 1:npcs, verbose = verbose)

    adj = obj.cluster@graphs$RNA_snn

    return(adj)
}

