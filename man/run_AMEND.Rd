% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_AMEND.R
\name{run_AMEND}
\alias{run_AMEND}
\title{Identify active modules from an interaction network}
\usage{
run_AMEND(
  graph = NULL,
  adj_matrix = NULL,
  edge_list = NULL,
  n = 50,
  data = NULL,
  node_type = NULL,
  brw.attr = NULL,
  FUN = NULL,
  FUN.params = NULL,
  heterogeneous = FALSE,
  multiplex = FALSE,
  aggregate.multiplex = NULL,
  normalize = c("degree", "modified_degree"),
  k = 0.5,
  degree.bias = NULL,
  jump.prob = NULL,
  net.weight = NULL,
  switch.layer.prob = NULL,
  layer.weight = NULL,
  verbose = FALSE,
  eta = NULL,
  identifier = 1,
  in.parallel = FALSE,
  n.cores
)
}
\arguments{
\item{graph, adj_matrix, edge_list}{A single graph like object in the form of an igraph, adjacency matrix, or edge list. Or, a named list containing multiple graph-like objects of the same type, to be merged. The merged graph must be connected. If a third column is provided in an edge list, these are taken as edge weights. Only one of graph, adj_matrix, or edge_list should be given, with priority given to graph, then adj_matrix, then edge_list. See 'Details' for correct naming conventions.}

\item{n}{Size of the final module to be approximated.}

\item{data}{A named list of named numeric vectors (list elements correspond to graph components), a named numeric vector, or a character scalar (denoting a vertex attribute of the input igraph object) containing the experimental data from which seed values for RWR will be derived according to \emph{FUN} and \emph{FUN.params} args. See 'Details' for correct naming conventions.}

\item{node_type}{A named list of named character vectors (list elements correspond to graph components), a named character vector, a character scalar (denoting a vertex attribute of the input igraph object), or NULL. Denotes the component and the multiplex layer of each node. If NULL and multiplex and/or heterogeneous, node labels must follow 'name|type_layer' naming scheme (e.g., MYC|gene_1). See 'Details' for correct naming conventions.}

\item{brw.attr}{A named list of named numeric vectors (list elements correspond to graph components), a named numeric vector, a character vector (of nodes of interest), a character scalar (denoting a vertex attribute of the input igraph object), or NULL. Biased random walk vertex attribute values should be non-negative, with larger values increasing the transition probabilities to a node in RWR. If NULL, all nodes are given a value of 1. Character vector gives nodes of interest towards which transition probabilities are increased. See 'Details' for biased random walk info.}

\item{FUN}{A function, named list of functions, named list of character scalars, a single character scalar, or NULL. Function for transforming values in \emph{data} to derive seed values for RWR. Names correspond to graph components. Character strings should correspond to a default function: one of 'binary', 'shift_scale', 'p_value', or 'exp'. NULL means no transformation is done to values in \emph{data}. See 'Details' for descriptions of default functions.}

\item{FUN.params}{A named list of lists of named function arguments, a named list of named function arguments, or NULL. Function arguments to be passed to \emph{FUN}. Names should match names in \emph{FUN}.}

\item{heterogeneous}{Logical. If TRUE, graph is considered heterogeneous (more than one distinct node type, e.g., proteins and metabolites), and \emph{node_type} must be included as an argument or graph vertex attribute.}

\item{multiplex}{Logical. If true, graph is assumed to contain multiplex components.}

\item{aggregate.multiplex}{A named list (for aggregating) or NULL (for no aggregating). The list element 'primary' contains the name of the primary layer for a multiplex component whose edges will be used during subnetwork identification. The multiplex component is collapsed onto this primary layer. The list element 'agg.method' contains a character scalar referring to an aggregation function for aggregating vertex attributes ('mean', 'median', 'sum', 'gmean', 'hmean'). gmean=geometric, hmean=harmonic.}

\item{normalize}{Normalization scheme of adjacency matrix for random walk with restart}

\item{k}{Value between 0 and 1. When normalize = "modified_degree", the adjacency matrix is first left and right multiplied by a diagonal matrix of node degrees, which is raised to the power -\emph{k}. As \emph{k} increases, edge weights are penalized more for the degrees of their adjacent nodes.}

\item{degree.bias}{List or NULL. NULL (default) for no degree bias adjustment. List names should include 'component', which is a character string of length greater than or equal to 1 specifying the components or layers to which the degree bias adjustment will be applied, and 'method', which is a character string of 'BS' (bistochastic scaling), 'IN' (inflation-normalization), or 'SDS' (stationary distribution scaling).}

\item{jump.prob, net.weight}{A named vector, or NULL. \emph{jump.prob} is the probability of random walker jumping from one component of graph to another in RWR. \emph{net.weight} is the relative weight given to nodes of a component of graph, applied to seed vector in RWR. Only used when heterogeneous=TRUE.}

\item{switch.layer.prob, layer.weight}{A named list of named vectors, or NULL. \emph{switch.layer.prob} is the probability of random walker to switch from current layer in a multiplex to another layer in same component. \emph{layer.weight} is relative weight given to nodes of a layer of a component of graph, applied to seed vector in RWR. List element names correspond to multiplex components, and vector names correspond to layers within a multiplex.}

\item{verbose}{Logical. Whether to output current iteration number to show progress.}

\item{eta}{Starting filtering rate. If NULL (default), a value is chosen based on the input network size and parameter \emph{n}.}

\item{identifier}{Optional. For use when performing many runs of AMEND to keep track of progress.}

\item{in.parallel}{Logical. Run computations in parallel.}

\item{n.cores}{number of cores to use in parallel computations}
}
\value{
a named list with the following elements:
\itemize{
\item module: the final module (i.e., subnetwork)
\item score: final module score
\item subnetworks: a list of node names contained in intermediate subnetworks
\item stats: network statistics
\item time: run time
\item input_params: list of input parameters
}
}
\description{
Identifies an active module, i.e., subnetwork with large experimental values, through an iterative optimization procedure using RWR, which can accommodate multiplex/heterogeneous networks.
}
\details{
Given an interaction network and experimental scores, AMEND attempts to find an active module by iteratively applying random walk with restart (RWR) and a heurstic solution to a maximum-weight connected subgraph problem.

Briefly, RWR is applied to the network to get diffusion scores which are shifted by some percentile (filtering rate/eta, which decreases exponentially each iteration). The network with diffusion scores are input into \code{heinz()}, an algorithm that heuristically finds a subnetwork with maximum node weight. This subnetwork is assigned a score based on its mean experimental values and mean core-clustering coefficients and then input into the next iteration. The algorithm stops when there is no change between iterations, and the subnetwork with maximum score is returned.

This function can also accommodate multiplex and/or heterogeneous networks, i.e., networks with more than one node type (e.g., proteins and metabolites) and more than one edge/data type within a node type (e.g., functional vs. physical interactions, or proteomic & transcriptomic data). This information is incorporated in RWR for multiplex/heterogeneous graphs (RWR-MH). This introduces the \emph{jump.prob}, \emph{net.weight}, \emph{switch.layer.prob}, and \emph{layer.weight} parameters. Importantly, this allows for seamless multi-omic integration in the context of active module identification.

See \code{create_integrated_graph()} details for proper naming conventions for node labels and function arguments when input network is multiplex and/or heterogeneous.

\emph{FUN} and \emph{FUN.params} are used to calculate seed values (which must be non-negative) from the values in \emph{data} for use in \code{RWR()}. There are default functions you can specify using character scalars. For NULL, no transformation is done. 'binary' assumes that all values in \emph{data} are 0 or 1 so no transformation is done. 'shift_scale' is for data types with a range centered about zero and takes two parameters: \emph{DOI} (direction of interest: 1 for positive, -1 for negative, and 0 for both) and \emph{w} (numeric value between 0 and 1). It takes the absolute values of \emph{data} and then down-weights nodes that weren't in DOI by \emph{w} (if \emph{DOI}=0, \emph{w} is coerced to 1). 'p_value' assumes \emph{data} are p-values and calculates \emph{-log10()} . 'exp' assumes \emph{data} are log fold changes and has a \emph{DOI} arg. For \emph{DOI}=-1,1, it exponentiates the product of \emph{DOI} and \emph{logFC}. For \emph{DOI}=0, it exponentiates abs( \emph{logFC}).

\code{run_AMEND()} can also implement a biased random walk based on a user-defined vertex attribute \emph{brw.attr}. These should be non-negative with larger values increasing transition probabilities to a node in RWR. This allows you to integrate other information (e.g., survival data) in active module identification. Another possible use is to artificially center modules around nodes of interest by using a character vector for \emph{brw.attr}. This will assign node values that decrease as a function of distance from these nodes of interest. The biased random walk is done by left multiplying the adjacency matrix by a diagonal matrix whose elements are values given by \emph{brw.attr}.

The \emph{aggregate.multiplex} argument specifies whether multiplex layers should be collapsed for subnetwork identification. This involves running RWR on the full network and then collapsing the multiplex layers onto one layer given by the "primary" argument in the input list. RWR scores are aggregated using the 'agg.method' argument in input list. The subnetwork that is returned will contain the collapsed multiplex component rather than nodes from individual layers. This can be done for none, some, or all of the multiplex components of the input graph.

The \emph{degree.bias} argument specifies whether degree bias should be adjusted for. Degree bias arises when certain molecules are studied more often then others, or if certain technologies systematically favor the identification of interactors of molecules with certain characteristics. There are 3 methods available: Inflation-normalization ('IN'), Bistochastic Scaling ('BS'), and Stationary Distribution Scaling ('SDS').

See the \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10324253/}{manuscript} for more details on the AMEND algorithm. There have been modifications to the algorithm as presented in the original paper.
}
\examples{
# Attach igraph library
library(igraph)

# Inspect the igraph object included in AMEND package.
# One can see it has vertex attributes name, symbol, and ECI
glut4_graph
head(V(glut4_graph)$ECI)

\dontrun{
# Use run_AMEND() with an igraph object with a vertex attribute matching data arg
subnet1 = run_AMEND(graph = glut4_graph, n = 25, data = "ECI", FUN = "shift_scale",
                    FUN.params = list(DOI = -1, w = 0.5), normalize = "degree")

# Use run_AMEND() with an adjacency matrix and a vector of node scores
subnet2 = run_AMEND(adj_matrix = glut4_adjM, data = eci_scores, FUN = "shift_scale",
                    FUN.params = list(DOI = -1, w = 0.5), normalize = "degree")

#=== Monoplex-heterogeneous network ===#
# Using a list of edge lists, with data as list input
# Create object for 'data' argument
data.list.input = vector("list", 2)
names(data.list.input) = c("prot", "meta")
for(i in seq_along(data.list.input)){
  uniq.names = extract_string(unique(c(edgelists[[i]][,1], edgelists[[i]][,2])), "\\\\|", 1)
  data.list.input[[i]] = rnorm(n = length(uniq.names), mean = 0, sd = 1)
  names(data.list.input[[i]]) = uniq.names
}

# Setting multiplex/heterogeneous parameters
net.weight = c(prot = 0.5, meta = 0.5)
jump.prob = c(prot = 0.7, meta = 0.3)

# FUN and FUN.params
FUN = list(prot = "shift_scale", meta = function(x, k) abs(x) * k)
FUN.params = list(prot = list(DOI = 1, w = 0.5), meta = list(k = 1.2))
subnet3 = run_AMEND(edge_list = edgelists, n = 50, data = data.list.input, node_type = NULL,
                    brw.attr = NULL, FUN = FUN, FUN.params = FUN.params, heterogeneous = TRUE,
                    multiplex = FALSE, normalize = "degree", jump.prob = jump.prob,
                    net.weight = net.weight, verbose = TRUE)

#=== Multiplex-heterogeneous network ===#
# List of graphs, data as vertex attribute
# Create data vertex attributes, called 'scores'
v.attr.name = "scores"
node.types = c("prot_1", "prot_2", "prot_3", "meta")
for(i in seq_along(node.types)){
  igraph::vertex_attr(list.of.graphs[[node.types[i]]],
                      v.attr.name) = runif(vcount(list.of.graphs[[node.types[i]]]))
}

# Setting multiplex/heterogeneous parameters
layer.weight = list(prot = rep(1/3, 3))
switch.layer.prob = list(prot = c(prot_1 = 0.5, prot_2 = 0.2, prot_3 = 0.8))
net.weight = c(prot = 0.5, meta = 0.5)
jump.prob = c(prot = 0.5, meta = 0.3)

subnet4 = run_AMEND(graph = list.of.graphs, n = 100, data = "scores", node_type = NULL,
                    brw.attr = NULL, FUN = "p_value", FUN.params = NULL,
                    heterogeneous = TRUE, multiplex = TRUE, normalize = "modified_degree",
                    k = 0.5, jump.prob = jump.prob, net.weight = net.weight,
                    switch.layer.prob = switch.layer.prob,
                    layer.weight = layer.weight, verbose = TRUE)

#=== Multiplex-homogeneous network ===#
# List of adjacency matrices, data as list input
# Create a list of binary vectors (e.g., significantly up- or down-regulated genes)
data.list.input = vector("list", 3)
names(data.list.input) = c("prot_1", "prot_2", "prot_3")
for(i in seq_along(data.list.input)){
  data.list.input[[i]] = sample(x = c(0,1), size = nrow(list.of.adjmats[[i]]),
                                replace = T, prob = c(0.9, 0.1)) # binary
  names(data.list.input[[i]]) = rownames(list.of.adjmats[[i]])
}

# Setting multiplex/heterogeneous parameters
layer.weight = list(prot = c(prot_1 = 0.4, prot_2 = 0.3, prot_3 = 0.3))
switch.layer.prob = list(prot = c(prot_1 = 0.4, prot_2 = 0.5, prot_3 = 0.5))

subnet5 = run_AMEND(adj_matrix = list.of.adjmats, n = 100, data = data.list.input,
                    node_type = NULL, brw.attr = NULL, FUN = "binary",
                    FUN.params = NULL, heterogeneous = FALSE, multiplex = TRUE,
                    normalize = "modified_degree", k = 0.5, layer.weight = layer.weight,
                    switch.layer.prob = switch.layer.prob, verbose = TRUE)
}

}
\seealso{
\code{\link[=create_integrated_graph]{create_integrated_graph()}}, \code{\link[=list_input_checks]{list_input_checks()}}, \code{\link[=transition_matrix]{transition_matrix()}}, \code{\link[=RWR]{RWR()}}, \code{\link[=heinz]{heinz()}}, \code{\link[=get_subnetwork]{get_subnetwork()}}, \code{\link[=inflate_normalize]{inflate_normalize()}}, \code{\link[=bistochastic_scaling]{bistochastic_scaling()}}
}
