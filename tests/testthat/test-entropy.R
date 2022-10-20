# This tests the perCellEntropy() function.
# library(testthat); library(TSCAN); source("test-entropy.R")

set.seed(1019201)

test_that("perCellEntropy works as expected", {
    x <- matrix(rpois(10000, lambda=5), ncol=200)
    ent <- perCellEntropy(x)
    expect_identical(length(ent), ncol(x))

    # Greater dispersion = less homogeneous = lower entropy.
    y <- matrix(rnbinom(10000, mu=5, size=0.5), ncol=200)
    ent2 <- perCellEntropy(y)
    expect_true(mean(ent) > mean(ent2))
})

test_that("perCellEntropy works with different matrix representations", {
    x <- matrix(rpois(10000, lambda=0.5), ncol=200)
    ent <- perCellEntropy(x)

    library(Matrix)
    y <- as(x, "dgCMatrix")
    ent2 <- perCellEntropy(y)
    expect_identical(ent, ent2)

    library(DelayedArray)
    z <- DelayedArray(x)
    ent3 <- perCellEntropy(z)
    expect_identical(ent, ent3)
})

test_that("perCellEntropy is robust to all-zero columns", {
    x <- matrix(rpois(10000, lambda=0.5), ncol=200)
    x[,1] <- 0

    ent <- perCellEntropy(x)
    expect_identical(length(ent), ncol(x))
    expect_identical(ent[1], NA_real_)
})
