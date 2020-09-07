# This tests the createClusterMST function.
# library(testthat); library(TSCAN); source("test-create.R")

set.seed(1000)

test_that("MST construction works as expected", {
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(0, 4)) 
    mst <- createClusterMST(y, cluster=NULL)

    vertices <- names(igraph::V(mst))
    expect_identical(vertices, rownames(y))
    expect_identical(vertices[igraph::degree(mst)==1], c("A", "D"))
    expect_identical(vertices[igraph::degree(mst)==2], c("B", "C"))

    # For a more complex example.
    y <- rbind(A=c(0, 1), B=c(0, 0), C=c(1, 1), D=c(-1, 1)) 
    mst <- createClusterMST(y, cluster=NULL)

    vertices <- names(igraph::V(mst))
    expect_identical(vertices, rownames(y))
    expect_identical(vertices[igraph::degree(mst)==1], c("B", "C", "D"))
    expect_identical(vertices[igraph::degree(mst)==3], "A")
})

test_that("MST edge and vertex attributes make sense", {
    y <- matrix(rnorm(20), ncol=2)
    rownames(y) <- LETTERS[1:10]
    mst <- createClusterMST(y, cluster=NULL)

    # Check that the distances between edges are what is expected.
    stuff <- Matrix::which(mst[] > 0, arr.ind=TRUE)
    coords <- igraph::V(mst)$coordinates

    for (i in seq_len(nrow(stuff))) {
        M <- stuff[i,1]
        N <- stuff[i,2]
        expect_equal(mst[M, N], sqrt(sum((coords[[M]] - coords[[N]])^2)))
    }

    # Centers must be named!
    y <- matrix(rnorm(20), ncol=2)
    expect_error(createClusterMST(y, cluster=NULL), "must have row names")
})

test_that("MST construction works as expected with clusters", {
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(0, 4))
    mst <- createClusterMST(y, cluster=NULL)

    clusters <- sample(rownames(y), 1000, replace=TRUE)
    y0 <- y[clusters,,drop=FALSE] + runif(1000, -0.1, 0.1)
    mst0 <- createClusterMST(y0, clusters=clusters)
    expect_identical(mst[] > 0, mst0[] > 0)

    # For a more complex example.
    y <- rbind(A=c(0, 1), B=c(0, 0), C=c(1, 1), D=c(-1, 1))
    mst <- createClusterMST(y, cluster=NULL)

    clusters <- sample(rownames(y), 1000, replace=TRUE)
    y0 <- y[clusters,,drop=FALSE] + runif(1000, -0.1, 0.1)
    mst0 <- createClusterMST(y0, clusters=clusters)
    expect_identical(mst[] > 0, mst0[] > 0)
})

test_that("MST construction works as expected with columns", {
    y <- matrix(rnorm(100), ncol=10)
    rownames(y) <- 1:10
    mst <- createClusterMST(y[,1:5], cluster=NULL)
    mst0 <- createClusterMST(y, cluster=NULL, columns=1:5)
    expect_identical(mst[], mst0[])
})

test_that("MST construction works as expected with outgroup specification", {
    # No effect with default outgroup settings.
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(0, 4)) 
    ref <- createClusterMST(y, cluster=NULL)
    mst <- createClusterMST(y, outgroup=TRUE, cluster=NULL)
    expect_identical(ref[], mst[])

    # Maximal effect with crazy outgroup settings.
    mst <- createClusterMST(y, outgroup=1e-8, cluster=NULL)
    expect_identical(dim(mst[]), c(4L, 4L))
    expect_identical(sum(mst[]), 0)

    # Default effects.
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(10, 3), D=c(10, 4))
    mst <- createClusterMST(y, outgroup=TRUE, cluster=NULL)
    expect_identical(igraph::components(mst)$no, 2L)
    expect_false(igraph::are_adjacent(mst, "B", "C"))

    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(10, 4), E=c(0,4)) 
    mst <- createClusterMST(y, outgroup=TRUE, cluster=NULL)
    expect_identical(sum(mst[]["D",]), 0)
    expect_identical(sum(mst[][,"D"]), 0)
})

library(SingleCellExperiment)
test_that("MST construction works with SE and SCE inputs", {
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(0, 4)) 
    ref <- createClusterMST(y, cluster=NULL)

    se <- SummarizedExperiment(t(y))
    mst <- createClusterMST(se, cluster=NULL, assay.type=1)
    expect_identical(ref[], mst[])

    sce <- SingleCellExperiment(t(y))
    mst <- createClusterMST(sce, cluster=NULL, assay.type=1)
    expect_identical(ref[], mst[])

    y2 <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(10, 4))
    reducedDim(sce, "PCA") <- y2
    mst <- createClusterMST(sce, cluster=NULL, use.dimred="PCA")
    expect_false(identical(ref[], mst[]))
})
