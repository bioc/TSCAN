# This tests the edge reporting function
# library(testthat); library(TSCAN); source("test-report.R")

test_that("MST segment reporting works as expected", {
    y <- rbind(A=c(0, 1), B=c(0, 2), C=c(0, 3), D=c(0, 4)) 
    mst <- createClusterMST(y, clusters=NULL)

    out <- reportEdges(y, mst, clusters=NULL)
    expect_identical(colnames(out), c("edge", "dim1", "dim2"))
    expect_identical(nrow(out), length(igraph::E(mst)) * 2L)
    
    alt <- reportEdges(x=NULL, mst) # same as re-using 'y'.
    expect_identical(out, alt)

    alt <- reportEdges(y[rev(rownames(y)),,drop=FALSE], mst, clusters=NULL) # robust to shuffling
    expect_identical(out, alt)

    named <- y
    colnames(named) <- c("PC1", "PC2")
    out2 <- reportEdges(named, mst, clusters=NULL)
    expect_identical(colnames(out2), c("edge", "PC1", "PC2"))

    out <- reportEdges(y, mst, combined=FALSE, clusters=NULL)
    expect_identical(colnames(out$start), c("edge", "dim1", "dim2"))
    expect_identical(colnames(out$end), c("edge", "dim1", "dim2"))
    expect_identical(nrow(out$start), length(igraph::E(mst)))
    expect_identical(nrow(out$end), length(igraph::E(mst)))
})

set.seed(100100)
test_that("MST segment reporting works with clusters= and columns=", {
    stuff <- matrix(rnorm(1000), ncol=10)
    clust <- kmeans(stuff,5)$cluster

    mst <- createClusterMST(stuff, clusters=clust)
    out <- reportEdges(stuff, mst, clusters=clust)
    expect_identical(nrow(out), length(igraph::E(mst)) * 2L)

    centered <- rowmean(stuff, clust)
    out2 <- reportEdges(centered, mst, clusters=NULL)
    expect_identical(out, out2)

    stuffx <- matrix(rnorm(1000), ncol=10)
    combined <- cbind(stuff, stuffx)
    out3 <- reportEdges(combined, mst, clusters=clust, columns=1:10)
    expect_identical(out, out3)
})

set.seed(100100)
test_that("MST segment reporting works with SE and SCE objects", {
    y <- matrix(rnorm(1000), ncol=10)
    clust <- kmeans(y,5)$cluster

    library(SingleCellExperiment)
    se <- SummarizedExperiment(t(y))
    mst <- createClusterMST(se, cluster=clust, assay.type=1)
    rep <- reportEdges(se, mst, cluster=clust, assay.type=1)
    expect_identical(rep, reportEdges(x=NULL, mst))

    sce <- SingleCellExperiment(t(y))
    rep2 <- reportEdges(sce, mst, cluster=clust, assay.type=1)
    expect_identical(rep, rep2)

    assay(sce) <- 2*assay(sce) # check it isn't just pulling it out of the assays.
    reducedDim(sce, "PCA") <- y
    rep3 <- reportEdges(sce, mst, cluster=clust, use.dimred="PCA")
    expect_identical(rep, rep3)
})

