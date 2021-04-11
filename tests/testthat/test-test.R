# This tests the testPseudotime functionality.
# library(testthat); library(TSCAN); source("test-test.R")

set.seed(01001)

test_that("basic tests work as expected", {
    y <- matrix(rnorm(10000), ncol=100)
    u <- runif(100)

    out <- testPseudotime(y, u)
    expect_type(out$logFC, "double")

    out <- testPseudotime(y, u, get.spline.coef=TRUE)
    expect_type(out$spline1, "double")

    test <- rbind(rnorm(100))
    x <- rnorm(100)
    ref <- testPseudotime(test, x)
    expect_true(ref$p.value > 0.1)

    test <- rbind(jitter(1:100))
    x <- jitter(1:100)
    ref <- testPseudotime(test, x)
    expect_true(ref$p.value < 0.01)
})

test_that("handles large numbers of duplicated pseudotime values", {
    y <- matrix(rnorm(10000), ncol=100)
    u <- c(numeric(50), runif(50))
    expect_error(out <- testPseudotime(y, u), NA)
    expect_error(out <- testPseudotime(y, u, df=51), "not enough unique")
})

test_that("handles NA pseudotime values", {
    y <- matrix(rnorm(10000), ncol=100)
    u <- c(numeric(50), runif(50))
    u[1:10] <- NA

    out <- testPseudotime(y, u)
    ref <- testPseudotime(y[,-(1:10)], u[-(1:10)])
    expect_identical(out, ref)
})

set.seed(01002)

test_that("tests handle blocking", {
    y <- matrix(rnorm(10000), ncol=100)
    u <- runif(100)
    b <- sample(3, ncol(y), replace=TRUE)

    out <- testPseudotime(y, u, block=b)

    expect_identical(length(out$per.block), 3L)

    a1 <- testPseudotime(y[,b==1], u[b==1])
    a2 <- testPseudotime(y[,b==2], u[b==2])
    a3 <- testPseudotime(y[,b==3], u[b==3])
    expect_equivalent(out$per.block[,1], a1)

    refp <- metapod::combineParallelPValues(list(a1$p.value, a2$p.value, a3$p.value), method="stouffer", weights=table(b))
    expect_identical(out$p.value, refp$p.value)

    # Works with parallelization.
    library(BiocParallel)
    alt <- testPseudotime(y, u, block=b, BPPARAM=SnowParam(2))
    expect_identical(out, alt)

    # With and without blocking works as expected.
    test <- rbind(jitter(b))
    x <- jitter(b)
    ref <- testPseudotime(test, x)
    expect_true(ref$p.value <= 0.01)
    out <- testPseudotime(test, x, block=b)
    expect_true(out$p.value > 0.1)
})

test_that("tests handle multiple pseudotime values", {
    y <- matrix(rnorm(10000), ncol=100)
    u <- matrix(runif(200), ncol=2)
    u[1:50,1] <- NA
    u[51:100,2] <- NA

    out <- testPseudotime(y, u)

    ref1 <- testPseudotime(y, u[,1])
    expect_identical(ref1, out[[1]])

    ref2 <- testPseudotime(y, u[,2])
    expect_identical(ref2, out[[2]])
})

test_that("testPseudotime supports addition of arbitrary row.data", {
    y <- matrix(rnorm(10000), ncol=100)
    u <- c(numeric(50), runif(50))
    u[1:10] <- NA

    rownames(y) <- paste0("GENE_", seq_len(nrow(y)))
    df <- DataFrame(Symbol=sample(letters, nrow(y), replace=TRUE))
    out <- testPseudotime(y, u, row.data=df)

    expect_identical(rownames(out), rownames(y))
    expect_identical(out$Symbol, df$Symbol)

    ref <- testPseudotime(y, u)
    expect_identical(ref, out[,-1])
})

