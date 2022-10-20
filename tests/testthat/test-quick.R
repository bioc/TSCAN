# This tests the quickPseudotime wrapper function.
# library(testthat); library(TSCAN); source("test-quick.R")

set.seed(1001010)
ncells <- 100
u <- matrix(rpois(20000, 5), ncol=ncells)
pca <- matrix(runif(ncells*5), ncells)
tsne <- matrix(rnorm(ncells*2), ncells)
clusters <- kmeans(pca, 3)$cluster

test_that("quickPseudotime works as expected", {
    qout <- quickPseudotime(pca, clusters=clusters, others=list(PCA=pca, TSNE=tsne))
    expect_identical(names(qout$centered), c("PCA", "TSNE"))

    # Options have an effect.
    qout2 <- quickPseudotime(pca, clusters=clusters, others=list(PCA=pca, TSNE=tsne), with.mnn=TRUE)
    expect_false(isTRUE(all.equal(qout$mst[], qout2$mst[])))

    ref <- quickPseudotime(pca[,1:2], clusters=clusters, others=list(PCA=pca, TSNE=tsne))
    out <- quickPseudotime(pca, clusters=clusters, others=list(PCA=pca, TSNE=tsne), columns=1:2)
    expect_identical(ref$ordering, out$ordering)
})

test_that("quickPseudotime works as expected for SE inputs", {
    ref <- quickPseudotime(pca, clusters=clusters, others=list(PCA=pca, TSNE=tsne))

    library(SingleCellExperiment)
    sce <- SingleCellExperiment(assays=list(counts=u), reducedDims=SimpleList(PCA=pca, tSNE=tsne))
    colLabels(sce) <- clusters
    out <- quickPseudotime(sce, use.dimred="PCA")

    expect_identical(ref$mst[], out$mst[])
    expect_identical(ref$ordering, out$ordering)

    # Also works for SE's.
    se <- SummarizedExperiment(t(pca))
    limited <- quickPseudotime(se, clusters=clusters, assay.type=1) 
    expect_identical(ref$mst[], limited$mst[])
    expect_identical(ref$ordering, limited$ordering)
})
