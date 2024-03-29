% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{list.of.graphs}
\alias{list.of.graphs}
\title{List of igraphs representing a multiplex-heterogeneous network}
\format{
A named list with 5 elements:
\itemize{
\item 'prot_1': igraph of layer 1 of protein multiplex component
\item 'prot_2': igraph of layer 2 of protein multiplex component
\item 'prot_3': igraph of layer 3 of protein multiplex component
\item 'meta': chemical-chemical interaction igraph
\item 'prot;meta': igraph of bipartite protein-chemical interaction network
}
}
\source{
Protein-protein interactions come from STRING. Chemical-chemical and chemical-protein interactions come from STITCH.
}
\usage{
list.of.graphs
}
\description{
5 separate igraph objects representing layers of a multiplex protein-protein interaction network, a chemical-chemical interaction network, and a bipartite chemical-protein interaction network. The multiplex layers were taken from 3 different clusters of the full PPI network after running the Louvain clustering algorithm.
}
\keyword{datasets}
