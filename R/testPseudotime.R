#' Test for differences along pseudotime
#'
#' Implements a simple method of testing for significant differences with respect to pseudotime,
#' based on fitting linear models with a spline basis matrix.
#'
#' @param x A numeric matrix-like object containing log-expression values for cells (columns) and genes (rows).
#' Alternatively, a \linkS4class{SummarizedExperiment} containing such a matrix.
#' @param pseudotime A numeric vector of length equal to the number of columns of \code{x}, containing the pseudotime orderings along a single lineage.
#' Alternatively, a single-column \linkS4class{PseudotimeOrdering} object containing such a vector in \code{\link{pathStat}(pseudotime)}.
#' @param df Integer scalar specifying the degrees of freedom for the splines.
#' @param get.lfc Logical scalar indicating whether to return an overall log-fold change along each path.
#' @param get.spline.coef Logical scalar indicating whether to return the estimates of the spline coefficients.
#' @param trend.only Logical scalar indicating whether only differences in the trend should be considered
#' when testing for differences between paths.
#' @param ... For the generic, further arguments to pass to specific method.
#' 
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' @param assay.type String or integer scalar specifying the assay containing the log-expression matrix.
#' @param block Factor of length equal to the number of cells in \code{x}, specifying the blocking factor.
#' @param BPPARAM A BiocParallelParam object from the \pkg{BiocParallel} package, used to control parallelization.
#'
#' @return
#' A \linkS4class{DataFrame} is returned containing the statistics for each gene (row),
#' including the p-value and its BH-adjusted equivalent.
#' If \code{get.lfc=TRUE}, an overall log-fold change is returned for each path.
#'
#' If \code{get.spline.coef=TRUE}, the estimated spline coefficients are also returned (single path)
#' or the differences in the spline fits to the first path are returned (multiple paths).
#' 
#' @details
#' Tis function fits a natural spline to the expression of each gene with respect to pseudotime.
#' It then does an ANOVA to test whether any of the spline coefficients are non-zero.
#' In this manner, genes exhibiting a significant (and potentially non-linear) trend
#' with respect to the pseudotime can be detected as those with low p-values.
#' 
#' For trajectories with multiple paths, only one path should be tested at a time.
#' This usually involves passing a single column of the matrix returned from \code{\link{orderCells}}.
#' Cells with \code{NA} values in \code{pseudotime} are assumed to be assigned to a different path and are ignored.
#'
#' By default, estimates of the spline coefficients are not returned as they are difficult to interpret.
#' Rather, a log-fold change of expression along each path is estimated
#' to provide some indication of the overall magnitude and direction of any change.
#'
#' \code{block} can be used to fit a separate linear model to each of multiple batches,
#' after which the statistics are combined across batches as described in \code{\link[scran]{testLinearModel}}.
#' This avoids potential confounding effects from batch-specific differences in the distribution of cells across pseudotime.
#' 
#' @author Aaron Lun
#'
#' @examples
#' y <- matrix(rnorm(10000), ncol=100)
#' u <- runif(100)
#' testPseudotime(y, u)
#'
#' # Handling a blocking factor.
#' b <- gl(2, 50)
#' testPseudotime(y, u, block=b)
#'
#' @seealso
#' \code{\link{orderCells}}, to generate the pseudotime matrix.
#'
#' \code{\link[scran]{testLinearModel}}, which performs the tests under the hood.
#'
#' @name testPseudotime
NULL

#' @importFrom stats p.adjust
#' @importFrom TrajectoryUtils pathStat
.test_pseudotime <- function(x, pseudotime, df=5, get.lfc=TRUE, get.spline.coef=FALSE, trend.only=TRUE, block=NULL, BPPARAM=NULL) {
    if (is(pseudotime, "PseudotimeOrdering")) {
        pseudotime <- pathStat(pseudotime)
    }

    output <- beachmat::rowBlockApply(x, FUN=.test_blocked_pseudotime, grid=TRUE, BPPARAM=BPPARAM, 
        pseudotime=pseudotime, df=df, get.lfc=get.lfc, get.spline.coef=get.spline.coef, block=block)

    output <- do.call(rbind, output)

    # Fixing FDRs due to row-level processing.
    output$FDR <- p.adjust(output$p.value, method="BH")
    for (i in seq_along(output$per.block)) {
        output$per.block[[i]]$FDR <- p.adjust(output$per.block[[i]]$p.value, method="BH")
    }

    output
}

#' @importFrom stats model.matrix
#' @importFrom S4Vectors metadata
.test_blocked_pseudotime <- function(x, pseudotime, ..., block) {
    pseudotime <- drop(pseudotime)

    if (is.null(block)) {
        .test_solo_pseudotime(x, pseudotime=pseudotime, ...)
    } else {
        by.block <- split(seq_len(ncol(x)), block)
        ncells <- lengths(by.block)
        valid <- logical(length(by.block))

        for (b in seq_along(by.block)) {
            current <- by.block[[b]]
            output <- .test_solo_pseudotime(x[,current,drop=FALSE], pseudotime=pseudotime[current], ...)
            by.block[[b]] <- output
            resid.df <- metadata(output)$residual.df
            valid[b] <- !is.na(resid.df) && resid.df > 0L
        }

        scran::combineBlocks(
            by.block,
            method="z", 
            geometric=FALSE,
            equiweight=FALSE, 
            weights=ncells,
            ave.fields=setdiff(colnames(by.block[[1]]), c("p.value", "FDR")),
            pval.field="p.value", 
            valid=valid
        )
    }
}

.test_solo_pseudotime <- function(x, pseudotime, df, get.lfc, get.spline.coef) {
    keep <- !is.na(pseudotime)
    pseudotime <- pseudotime[keep]

    x <- x[,keep,drop=FALSE] 
    design <- .forge_spline_basis_design(pseudotime, df=df)
    output <- scran::testLinearModel(x, design=design, coefs=2:ncol(design))
 
    if (get.lfc) {
        prior <- colnames(output)
        design.lfc <- model.matrix(~pseudotime)
        output$logFC <- scuttle::fitLinearModel(x, design=design.lfc)$coefficients[,2]
        output <- output[,c("logFC", prior)] 
    }

    if (!get.spline.coef) {
        output <- output[,setdiff(colnames(output), colnames(design))]
    }

    output
}

#' @importFrom stats model.matrix
.forge_spline_basis_design <- function(p, df) { 
    # Uniquify'ing to avoid non-full rank problems when
    # many of the quantiles are stacked at the same position. 
    up <- unique(p)
    if (length(up) <= df) {
        stop("'not enough unique pseudotime values for the specified 'df'")
    }

    basis <- splines::ns(up, df=df)
    colnames(basis) <- sprintf("spline%i", seq_len(df))
    basis <- basis[match(p, up),,drop=FALSE]

    cbind(Intercept=rep(1, length(p)), basis)
}

#' @export
#' @rdname testPseudotime
setGeneric("testPseudotime", function(x, ...) standardGeneric("testPseudotime"))

#' @export
#' @rdname testPseudotime
setMethod("testPseudotime", "ANY", .test_pseudotime)

#' @export
#' @rdname testPseudotime
#' @importFrom SummarizedExperiment assay
setMethod("testPseudotime", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .test_pseudotime(assay(x, assay.type), ...)
})
