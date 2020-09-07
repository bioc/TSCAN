# This tests the orderCells function.
# library(testthat); library(TSCAN); source("test-order.R")

set.seed(99000001)
test_that("MST ordering works as expected for straight lines", {
    centers <- rbind(A=c(1, 0), B=c(3, 0), C=c(4.5, 0), D=c(6, 0))
    clusters <- sample(rownames(centers), 1000, replace=TRUE)
    y <- centers[clusters,,drop=FALSE]
    y <- y + rnorm(length(y))

    mst <- createClusterMST(centers, clusters=NULL)
    mapping <- mapCellsToEdges(y, mst, clusters=clusters)
    ordering <- orderCells(mapping, mst)

    m <- match(clusters, rownames(centers))
    left.lim <- c(min(centers[,1]), centers[,1])[m]
    right.lim <- c(centers[,1], max(centers[,1]))[m+1]
    expect_equivalent(as.numeric(ordering), pmin(right.lim, pmax(left.lim, y[,1])) - 1)
})

set.seed(99000002)
test_that("MST ordering works as expected for branching", {
    centers <- rbind(base=c(0, 0), A=c(1, 0), B=c(0, 2), C=c(0, -3), D=c(-4, 0))
    clusters <- sample(rownames(centers), 1000, replace=TRUE)
    y <- centers[clusters,,drop=FALSE]
    y <- y + runif(length(y), -0.4, 0.4)

    mst <- createClusterMST(centers, clusters=NULL)
    mapping <- mapCellsToEdges(y, mst, clusters=clusters)
    ordering <- orderCells(mapping, mst, start="base")

    # Everyone gets assigned to only one pseudotime.
    expect_identical(ncol(ordering), nrow(centers) - 1L)
    expect_true(all(rowSums(!is.na(ordering))==1))

    # Computing the branched pseudo-times.
    expect_equivalent(ordering[clusters=="A",1], pmin(y[clusters=="A",1], 1))
    expect_equivalent(ordering[clusters=="B",2], pmin(y[clusters=="B",2], 2))
    expect_equivalent(ordering[clusters=="C",3], -pmax(y[clusters=="C",2], -3))
    expect_equivalent(ordering[clusters=="D",4], -pmax(y[clusters=="D",1], -4))
})

set.seed(99000002)
test_that("MST ordering works as expected for a complex case", {
    centers <- rbind(base=c(0, 0), A=c(1, 0), B=c(0, 2), C=c(0, -3), D=c(-4, 0))
    clusters <- sample(rownames(centers), 1000, replace=TRUE)
    y <- centers[clusters,,drop=FALSE]
    y <- y + runif(length(y), -0.4, 0.4)

    mst <- createClusterMST(centers, clusters=NULL)
    mapping <- mapCellsToEdges(y, mst, clusters=clusters)
    ordering <- orderCells(mapping, mst, start="base")

    # Everyone gets assigned to only one pseudotime.
    expect_identical(ncol(ordering), nrow(centers) - 1L)
    expect_true(all(rowSums(!is.na(ordering))==1))

    # Computing the branched pseudo-times.
    expect_equivalent(ordering[clusters=="A",1], pmin(y[clusters=="A",1], 1))
    expect_equivalent(ordering[clusters=="B",2], pmin(y[clusters=="B",2], 2))
    expect_equivalent(ordering[clusters=="C",3], -pmax(y[clusters=="C",2], -3))
    expect_equivalent(ordering[clusters=="D",4], -pmax(y[clusters=="D",1], -4))

    # Computing the expected partition for the base.
    base <- clusters=="base"
    expect_equivalent(
        rowMeans(ordering[base,,drop=FALSE], na.rm=TRUE),
        pmax(abs(y[base,1]), abs(y[base,2]))
    )
})

set.seed(990000021)
test_that("MST ordering shares values at zero", {
    X <- as.matrix(expand.grid(-1:1, -1:1))
    Y <- rbind(A=c(-1, -1), B=c(0, 0), C=c(1, -1))

    mst <- createClusterMST(Y, clusters=NULL)
    mapping <- mapCellsToEdges(X, mst, clusters=rep("B", nrow(X)))
    ordering <- orderCells(mapping, mst, start="B")

    expect_identical(sum(ordering==0, na.rm=TRUE)/ncol(ordering), 4)
    expect_identical(which(ordering[,1]==0), which(ordering[,2]==0))

    # Still the case when it's not happening on the starting node.
    Y2 <- rbind(Y, D=c(0, -1))
    mst <- igraph::make_graph(c("D", "B", "B", "A", "B", "C"), directed=FALSE)
    igraph::V(mst)$coordinates <- split(Y2, rep(rownames(Y2), ncol(Y2)))[names(igraph::V(mst))]
    igraph::E(mst)$weight <- c(1, sqrt(2), sqrt(2))

    mapping <- mapCellsToEdges(X, mst, clusters=rep("B", nrow(X)))
    ordering2 <- orderCells(mapping, mst, start="D")

    expect_identical(sum(ordering2==0, na.rm=TRUE)/ncol(ordering2), 1)
    expect_identical(sum(ordering2==1, na.rm=TRUE)/ncol(ordering2), 4)

    expect_identical(which(ordering2[,1]==1), which(ordering2[,2]==1))
    expect_identical(which(ordering2[,1]==0), which(ordering2[,2]==0))
})

set.seed(99000003)
test_that("MST ordering is robust to transformations", {
    centers <- matrix(rnorm(20, sd=3), ncol=2)
    rownames(centers) <- letters[1:10]

    clusters <- sample(rownames(centers), 1000, replace=TRUE)
    y <- centers[clusters,,drop=FALSE]
    y <- y + rnorm(length(y), 0.1)

    mst <- createClusterMST(centers, clusters=NULL)
    map <- mapCellsToEdges(y, clusters=clusters, mst)
    ref <- orderCells(map, mst)

    # For your viewing pleasure:
    # plot(y[,1], y[,2], col=topo.colors(21)[cut(rowMeans(ref, na.rm=TRUE), 21)])
    # stuff <- connectClusterMST(centers, mst, combined=FALSE)
    # segments(stuff$start$dim1, stuff$start$dim2, stuff$end$dim1, stuff$end$dim2, lwd=5)

    # Applying various transformation.
    mst2 <- createClusterMST(centers * 2, clusters=NULL)
    map2 <- mapCellsToEdges(y*2, clusters=clusters, mst2)
    out2 <- orderCells(map2, mst2)
    expect_equivalent(ref*2, out2)

    mstp1 <- createClusterMST(centers + 1, clusters=NULL)
    mapp1 <- mapCellsToEdges(y+1, clusters=clusters, mstp1)
    outp1 <- orderCells(mapp1, mstp1)
    expect_equivalent(ref, outp1)

    rotation <- matrix(c(cos(pi/4), sin(pi/4), -sin(pi/4), cos(pi/4)), ncol=2)
    mstr <- createClusterMST(centers %*% rotation, clusters=NULL)
    mapr <- mapCellsToEdges(y %*% rotation, clusters=clusters, mstr)
    outr <- orderCells(mapr, mstr)
    expect_equivalent(ref, outr)
})

set.seed(99000004)
test_that("MST ordering behaves with multiple components", {
    centers <- matrix(rnorm(20, sd=3), ncol=2)
    rownames(centers) <- letters[1:10]

    clusters <- sample(rownames(centers), 1000, replace=TRUE)
    y <- centers[clusters,,drop=FALSE]
    y <- y + rnorm(length(y), 0.1)

    mst <- createClusterMST(centers, clusters=NULL)
    map <- mapCellsToEdges(y, mst, clusters=clusters)
    ref <- orderCells(map, mst)

    # Doubling everything and seeing if the system recovers the right pseudotimes for each component.
    doublement <- rbind(centers, centers+100)
    rownames(doublement) <- c(rownames(centers), paste0(rownames(centers), "0"))
    y2 <- rbind(y, y+100)
    clusters2 <- c(clusters, paste0(clusters, "0"))

    # Bumping up 'outscale' to account for fluctuations in 'centers'.
    mst2 <- createClusterMST(doublement, clusters=NULL, outgroup=TRUE, outscale=10)
    comp <- igraph::components(mst2)
    expect_identical(comp$no, 2L)

    map2 <- mapCellsToEdges(y2, mst2, clusters=clusters2)
    out <- orderCells(map2, mst2)

    original <- !grepl("0", colnames(out))
    first <- seq_along(clusters)
    second <- first + length(clusters)

    expect_equivalent(out[first,original], ref)
    expect_equivalent(out[first,original], out[second,!original])
    expect_true(all(is.na(out[second,original])))
    expect_true(all(is.na(out[first,!original])))

    # Test passing of 'start' values. 
    expect_error(orderCells(map2, mst2, start=c("A")), "must have one cluster")
    by.comp <- split(names(comp$membership), comp$membership)
    test.starts <- vapply(by.comp, head, n=1, "")
    expect_error(orderCells(map2, mst2, start=test.starts), NA)
})

set.seed(99000005)
test_that("MST ordering behaves with multiple components where one has no edges", {
    centers <- matrix(rep(1:4, 2), ncol=2)
    centers[4,] <- centers[4,] * 10
    rownames(centers) <- letters[1:4]

    clusters <- sample(rownames(centers), 1000, replace=TRUE)
    y <- centers[clusters,,drop=FALSE]
    y <- y + rnorm(length(y), 0.1)

    mst <- createClusterMST(centers, clusters=NULL, outgroup=TRUE)
    expect_equivalent(igraph::degree(mst, "d"), 0L)

    map <- mapCellsToEdges(y, mst, clusters=clusters)
    ref <- orderCells(map, mst)

    expect_identical(ncol(ref), 2L)    
    expect_identical(sum(ref[,"d"], na.rm=TRUE), 0)
})
