#' Quick MST-based pseudotime
#'
#' A convenience wrapper to quickly compute a minimum spanning tree (MST) on the cluster centroids
#' to obtain a pseudotime ordering of the cells.
#'
#' @param x A numeric matrix of coordinates where each row represents a cell/sample and each column represents a dimension
#' (usually a PC or another low-dimensional embedding, but features or genes can also be used).
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} object
#' containing such a matrix in its \code{\link{assays}}, as specified by \code{assay.type}.
#' This will be transposed prior to use.
#'
#' Alternatively, for \linkS4class{SingleCellExperiment}s, this matrix may be extracted from its \code{\link{reducedDims}},
#' based on the \code{use.dimred} specification.
#' In this case, no transposition is performed.
#' @param clusters A vector or factor of length equal to the number of cells in \code{x},
#' specifying the cluster assignment for each cell.
#' @param start Passed to \code{\link{orderCells}}.
#' @param others List of numeric matrices with the same number of rows as \code{x}, to be passed to \code{\link{reportEdges}}.
#' This typically contains dimensionality reduction results, for use in visualizing the edges of the MST.
#' If \code{NULL}, defaults to a list containing \code{x}.
#' @param other.dimreds Logical scalar indicating whether all dimensionality reduction results in \code{x} 
#' should be appended onto the \code{others} list.
#' @param assay.type An integer or string specifying the assay to use from a SummarizedExperiment \code{x}.
#' @param outscale Passed to \code{\link{createClusterMST}}.
#' A different default is used here for historical reasons.
#' @param ... For the generic, further arguments to pass to the specific methods.
#'
#' For the ANY method, further arguments to pass to \code{\link{createClusterMST}}.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method
#' (if \code{use.dimred} is specified) or the ANY method (otherwise).
#'
#' @details
#' This function simply calls, in order:
#' \itemize{
#' \item \code{\link{rowmean}}, to compute the average low-dimensional coordinates for each cluster.
#' \item \code{\link{createClusterMST}} on the average coordinates created from \code{x}.
#' \item \code{\link{reportEdges}} on the average coordinates for all entries of \code{other}.
#' \item \code{\link{mapCellsToEdges}} on the per-cell coordinates in \code{x} with the constructed MST.
#' \item \code{\link{orderCells}} on the mappings generated from \code{x} onto the MST.
#' }
#'
#' @return
#' A \linkS4class{List} containing:
#' \itemize{
#' \item \code{centered}, a list of numeric matrices containing the averaged coordinates for each cluster.
#' Each matrix corresponds to a dimensionality reduction result in \code{x}.
#' \item \code{mst}, a \link{graph} object containing the cluster-level MST computed on the coordinates from \code{use}.
#' \item \code{ordering}, a \linkS4class{PseudotimeOrdering} object containing the ordering for various paths through the MST computed from \code{use}.
#' \item \code{connected}, a list of data.frames containing the edge coordinates between centers.
#' Each data.frame corresponds to a dimensionality reduction result in \code{x}.
#' }
#'
#' @seealso
#' \code{\link{createClusterMST}} and friends, for the functions that do the actual work.
#'
#' @author Aaron Lun
#' @examples
#' # Mocking up an SCE object:
#' ncells <- 100
#' u <- matrix(rpois(20000, 5), ncol=ncells)
#' pca <- matrix(runif(ncells*5), ncells)
#' tsne <- matrix(rnorm(ncells*2), ncells)
#'
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays=list(counts=u),
#'     reducedDims=SimpleList(PCA=pca, tSNE=tsne))
#'
#' # Clustering on our pretend PCA values:
#' clusters <- kmeans(pca, 3)$cluster
#'
#' # Quickly computing the pseudotime:
#' out <- quickPseudotime(sce, clusters, use.dimred="PCA")
#' out$mst
#' head(out$ordering)
#'
#' @name quickPseudotime
NULL

#' @importFrom S4Vectors List
#' @importFrom igraph V
.quick_pseudotime <- function(x, clusters, others=NULL, ..., outscale=1.5, start=NULL, columns=NULL) {
    # Do not center beforehand, as people might want to use with.mnn=TRUE.
    mst <- createClusterMST(x, clusters=clusters, ..., outscale=outscale, columns=columns)

    to.use <- lapply(others, rowmean, group=clusters)
    connected <- lapply(to.use, FUN=reportEdges, clusters=NULL, mst=mst, columns=columns)

    mapping <- mapCellsToEdges(x, clusters=clusters, mst=mst, columns=columns)
    ordering <- orderCells(mapping, mst, start=start)

    List(
        centered=to.use,
        mst=mst,
        ordering=ordering,
        connected=connected
    )
}

#' @export
#' @rdname quickPseudotime
setGeneric("quickPseudotime", function(x, ...) standardGeneric("quickPseudotime"))

#' @export
#' @rdname quickPseudotime
setMethod("quickPseudotime", "ANY", .quick_pseudotime)

#' @export
#' @rdname quickPseudotime
setMethod("quickPseudotime", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .quick_pseudotime(t(assay(x, assay.type)), ...)
})

#' @export
#' @rdname quickPseudotime
#' @importFrom SingleCellExperiment colLabels reducedDims reducedDim
setMethod("quickPseudotime", "SingleCellExperiment", function(x, clusters=colLabels(x, onAbsence="error"), 
    ..., others=NULL, use.dimred=NULL, other.dimreds=TRUE) 
{
    if (other.dimreds) {
        others <- c(others, as.list(reducedDims(x)))
    }
    if (!is.null(use.dimred)) {
        .quick_pseudotime(reducedDim(x, use.dimred), clusters=clusters, others=others, ...)
    } else {
        callNextMethod(x, clusters=clusters, others=others, ...)
    }
})
