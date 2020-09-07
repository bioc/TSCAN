# This tests the mapCellsToEdges function.
# library(testthat); library(TSCAN); source('test-mapping.R')

set.seed(10001)
test_that("mapCellsToEdges works correctly in general", {
    y <- matrix(rnorm(1000), ncol=10)
    clust <- kmeans(y,5)$cluster

    mst <- createClusterMST(y, cluster=clust)
    mapping <- mapCellsToEdges(y, mst, clusters=clust)
    expect_true(all(mapping$left.cluster == clust | mapping$right.cluster ==clust))

    # Expect the edge lengths to match up.
    mat <- mst[]
    edge.lens <- mat[cbind(
        match(mapping$left.cluster, colnames(mat)),
        match(mapping$right.cluster, colnames(mat))
    )]
    expect_equal(edge.lens, mapping$left.distance + mapping$right.distance)

    # Expect distances along the MST to be shorter than the actual distances.
    centered <- rowmean(y, clust)
    expect_true(all(
        mapping$left.distance <= sqrt(rowSums((y - centered[mapping$left.cluster,])^2))
    ))
    expect_true(all(
        mapping$right.distance <= sqrt(rowSums((y - centered[mapping$right.cluster,])^2))
    ))
})

test_that("mapCellsToEdges maps elements onto vertices correctly", {
    y <- matrix(rnorm(100), ncol=10)
    rownames(y) <- 1:10
    mst <- createClusterMST(y, cluster=NULL)
    mapping <- mapCellsToEdges(y, mst, clusters=rownames(y))

    expect_true(all(
        (mapping$left.cluster==1:10 & mapping$left.distance < 1e-8) |
        (mapping$right.cluster==1:10 & mapping$right.distance < 1e-8)
    ))

    # Works on extremes correctly.
    y <- matrix(1:10, nrow=10, ncol=2)
    rownames(y) <- 1:10
    mst <- createClusterMST(y, cluster=NULL)
    mapping <- mapCellsToEdges(rbind(c(0,0), c(100, 100)), mst, clusters=NULL)

    expect_true(mapping$right.cluster[1]==1 && mapping$right.distance[1]==0)
    expect_true(mapping$left.cluster[2]==10 && mapping$left.distance[2]==0)
})

set.seed(10002)
test_that("mapCellsToEdges handles free mapping correctly", {
    y <- matrix(rnorm(1000), ncol=10)
    clust <- kmeans(y,5)$cluster
    mst <- createClusterMST(y, cluster=clust)

    y2 <- matrix(rnorm(1000), ncol=10)
    free <- mapCellsToEdges(y2, mst, clusters=NULL)

    # Expect the edge lengths to match up.
    mat <- mst[]
    edge.lens <- mat[cbind(
        match(free$left.cluster, colnames(mat)),
        match(free$right.cluster, colnames(mat))
    )]
    expect_equal(edge.lens, free$left.distance + free$right.distance)

    # Expect the same behavior if we forced it to a cluster.
    forced1 <- mapCellsToEdges(y2, mst, clusters=free$left.cluster)
    expect_identical(free, forced1)

    forced2 <- mapCellsToEdges(y2, mst, clusters=free$right.cluster)
    expect_identical(free, forced2)
})

set.seed(100100)
test_that("mapCellsToEdges works with SE and SCE objects", {
    y <- matrix(rnorm(1000), ncol=10)
    clust <- kmeans(y,5)$cluster

    library(SingleCellExperiment)
    se <- SummarizedExperiment(t(y))
    mst <- createClusterMST(se, cluster=clust, assay.type=1)
    map <- mapCellsToEdges(se, mst, cluster=clust, assay.type=1)
    expect_identical(map, mapCellsToEdges(y, mst, cluster=clust))

    sce <- SingleCellExperiment(t(y))
    map2 <- mapCellsToEdges(sce, mst, cluster=clust, assay.type=1)
    expect_identical(map, map2)

    assay(sce) <- 2*assay(sce) # check it isn't just pulling it out of the assays.
    reducedDim(sce, "PCA") <- y
    map3 <- mapCellsToEdges(sce, mst, cluster=clust, use.dimred="PCA")
    expect_identical(map, map3)
})

