% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_AMEND.R
\name{set_brw_attr}
\alias{set_brw_attr}
\title{Set the Biased Random Walk attribute}
\usage{
set_brw_attr(graph, brw.attr, brw.attr.nm, brw.flag)
}
\arguments{
\item{graph}{Input graph}

\item{brw.attr}{Character vector. Denotes the nodes of interest (NOI) used to artificially seed the active module.}

\item{brw.attr.nm}{Character scalar. The name of the vertex attribute to be used for biased random walk.}

\item{brw.flag}{Logical. True returns graph with new biased random walk attribute. False returns input graph.}
}
\value{
igraph object
}
\description{
Calculates values for a biased random walk. This value is a decreasing function of minimum distance in the full graph to a node of interest (NOI): f(d,k) = exp(-k*d), where k >= 0, d=distance.
}
