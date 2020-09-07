\name{TSCANNews}
\title{TSCAN News}
\encoding{UTF-8}

\section{Version 1.28.0}{\itemize{
\item Added \code{createClusterMST()} to create a cluster-based MST from a variety of inputs,
migrated from the \pkg{scran} package.

\item Added \code{reportEdges()} to report edge coordinates for plotting.

\item Added \code{mapCellsToEdges()} to map cells to the closest edge on the MST.

\item Added \code{orderCells()} to compute a pseudotemporal ordering from mapped cells.

\item Added \code{quickPseudotime()} to wrap MST construction and ordering into a single call.

\item Added \code{testPseudotime()} to test for DE genes along one or more paths through a MST.

\item Added the \code{rowmean()} utility to compute column means for row groupings.
}}