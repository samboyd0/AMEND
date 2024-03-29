% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_AMEND.R
\name{edge_density_weighted}
\alias{edge_density_weighted}
\title{Calculate the edge density of a weighted graph}
\usage{
edge_density_weighted(graph, weight)
}
\arguments{
\item{graph}{Input graph}

\item{weight}{Logical. If true, the edge density is calculated taking into account edge weights}
}
\value{
weighted edge density of graph
}
\description{
\code{edge_density_weighted()} calculates the edge density of a graph, taking into account edge weights.
}
\details{
If weight is true, the edge density equals the sum of the edge weights divided by the maximum number of edges possible for the given graph. If weight is false, the edge density is the number of edges divided by the maximum number of edges possible for the given graph.
}
\examples{
AMEND:::edge_density_weighted(graph = glut4_graph, weight = TRUE)

AMEND:::edge_density_weighted(graph = glut4_graph, weight = FALSE)

# Compare to igraph::edge_density
igraph::edge_density(glut4_graph)

# Remove weight attribute
g = igraph::delete_edge_attr(glut4_graph, "weight")
AMEND:::edge_density_weighted(graph = g, weight = TRUE)

}
