#' Compute the per-cell entropy
#'
#' Compute the entropy of each cell, using this as a proxy for the differentiation status. 
#'
#' @param x A numeric matrix-like object containing counts for each cell (column) and feature (row).
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param BPPARAM A BiocParallelParam object from \pkg{BiocParallel}, specifying how calculations should be parallelized.
#' @param assay.type An integer or string specifying the assay to use from a SummarizedExperiment \code{x}.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' @details
#' Entropy values are computed from the proportion of counts assigned to each feature within a given cell.
#' The central idea is that undifferentiated cells have higher entropies because they are not yet committed to a single lineage,
#' and thus have low but persistent activity of the transcriptional programs for all lineages.
#' The cluster with the highest entropy values can be used to determine the \code{start} cluster in \code{\link{orderCells}}.
#' 
#' @return A numeric vector of entropies for all cells in \code{x}.
#' Cells with all-zero values in \code{x} will be assigned \code{NA} entropies.
#'
#' @examples
#' sce <- scuttle::mockSCE()
#' ent <- perCellEntropy(sce)
#' summary(ent)
#'
#' # Compute average entropy over mock clusters.
#' clusters <- sample(ncol(sce), 5)
#' by.cluster <- split(ent, clusters)
#' mean.cluster.ent <- vapply(by.cluster, mean, 0)
#'
#' @references
#' Grun D et al. (2016). 
#' De novo prediction of stem cell identity using single-cell transcriptome data.
#' \emph{Cell Stem Cell} 19, 266-77
#'
#' Gulati GS et al. (2020).
#' Single-cell transcriptional diversity is a hallmark of developmental potential.
#' \emph{Science} 367, 405-11
#'
#' Guo M et al. (2017)
#' SLICE: determining cell differentiation and lineage based on single cell entropy.
#' \emph{Nucleic Acids Res.} 45, e54
#'
#' @author Aaron Lun
#' @name perCellEntropy
NULL

#' @importFrom DelayedArray blockApply colAutoGrid
.per_cell_entropy <- function(x, BPPARAM=NULL) {
    output <- blockApply(x, FUN=.entropic_loop, as.sparse=NA, grid=colAutoGrid(x), BPPARAM=BPPARAM)
    unlist(output)
}

#' @importFrom Matrix colSums
#' @importClassesFrom DelayedArray SparseArraySeed
#' @importClassesFrom Matrix dgCMatrix
.entropic_loop <- function(x) {
    if (is(x, "SparseArraySeed")) {
        x <- as(x, "dgCMatrix")
    } 

    keep <- x > 0
    idx <- which(keep, arr.ind=TRUE)
    vals <- x[idx]

    totals <- colSums(x)
    p <- vals / totals[idx[,2]]
    entropy <- - p * log(p)

    f <- factor(idx[,2], levels=seq_len(ncol(x)))
    cell.entropy <- by(entropy, INDICES=f, FUN=sum)
    cell.entropy <- as.numeric(cell.entropy)

    cell.entropy[totals==0] <- NA_real_
    cell.entropy
}

#' @export
#' @rdname perCellEntropy
setGeneric("perCellEntropy", function(x, ...) standardGeneric("perCellEntropy"))

#' @export
#' @rdname perCellEntropy
setMethod("perCellEntropy", "ANY", .per_cell_entropy)

#' @export
#' @rdname perCellEntropy
setMethod("perCellEntropy", "SummarizedExperiment", function(x, ...,  assay.type="counts") {
    .per_cell_entropy(assay(x, assay.type), ...)
})
