#=======#
# To-Do:
#   - Create tests using testthat framework from R Packages by Hadley Wickham (12/27/23)
#   - Update functions documentation to describe the process of multiplex layer aggregation and the new biased random walk/nodes of interest functionality (1/29/24)
#=======#

#' @importFrom igraph vcount V V<- E E<- vertex_attr vertex_attr<-
#' @importFrom stats sd quantile
#' @importFrom foreach %dopar%

#' @title Identify active modules from an interaction network
#'
#' @description
#' Identifies an active module, i.e., subnetwork with large experimental values, through an iterative optimization procedure using RWR, which can accommodate multiplex/heterogeneous networks.
#'
#' @details
#'
#' Given an interaction network and experimental scores, AMEND attempts to find an active module by iteratively applying random walk with restart (RWR) and a heurstic solution to a maximum-weight connected subgraph problem.
#'
#' Briefly, RWR is applied to the network to get diffusion scores which are shifted by some percentile (filtering rate/eta, which decreases exponentially each iteration). The network with diffusion scores are input into `heinz()`, an algorithm that heuristically finds a subnetwork with maximum node weight. This subnetwork is assigned a score based on its mean experimental values and mean core-clustering coefficients and then input into the next iteration. The algorithm stops when there is no change between iterations, and the subnetwork with maximum score is returned.
#'
#' This function can also accommodate multiplex and/or heterogeneous networks, i.e., networks with more than one node type (e.g., proteins and metabolites) and more than one edge/data type within a node type (e.g., functional vs. physical interactions, or proteomic & transcriptomic data). This information is incorporated in RWR for multiplex/heterogeneous graphs (RWR-MH). This introduces the _jump.prob_, _net.weight_, _switch.layer.prob_, and _layer.weight_ parameters. Importantly, this allows for seamless multi-omic integration in the context of active module identification.
#'
#' See `create_integrated_graph()` details for proper naming conventions for node labels and function arguments when input network is multiplex and/or heterogeneous.
#'
#' _FUN_ and _FUN.params_ are used to calculate seed values (which must be non-negative) from the values in _data_ for use in `RWR()`. There are default functions you can specify using character scalars. For NULL, no transformation is done. 'binary' assumes that all values in _data_ are 0 or 1 so no transformation is done. 'shift_scale' is for data types with a range centered about zero and takes two parameters: _DOI_ (direction of interest: 1 for positive, -1 for negative, and 0 for both) and _w_ (numeric value between 0 and 1). It takes the absolute values of _data_ and then down-weights nodes that weren't in DOI by _w_ (if _DOI_=0, _w_ is coerced to 1). 'p_value' assumes _data_ are p-values and calculates _-log10()_ . 'exp' assumes _data_ are log fold changes and has a _DOI_ arg. For _DOI_=-1,1, it exponentiates the product of _DOI_ and _logFC_. For _DOI_=0, it exponentiates abs( _logFC_).
#'
#' `run_AMEND()` can also implement a biased random walk based on a user-defined vertex attribute _brw.attr_. These should be non-negative with larger values increasing transition probabilities to a node in RWR. This allows you to integrate other information (e.g., survival data) in active module identification. Another possible use is to artificially center modules around nodes of interest by using a character vector for _brw.attr_. This will assign node values that decrease as a function of distance from these nodes of interest. The biased random walk is done by left multiplying the adjacency matrix by a diagonal matrix whose elements are values given by _brw.attr_.
#'
#' The _aggregate.multiplex_ argument specifies whether multiplex layers should be collapsed for subnetwork identification. This involves running RWR on the full network and then collapsing the multiplex layers onto one layer given by the "primary" argument in the input list. RWR scores are aggregated using the 'agg.method' argument in input list. The subnetwork that is returned will contain the collapsed multiplex component rather than nodes from individual layers. This can be done for none, some, or all of the multiplex components of the input graph.
#'
#' The _degree.bias_ argument specifies whether degree bias should be adjusted for. Degree bias arises when certain molecules are studied more often then others, or if certain technologies systematically favor the identification of interactors of molecules with certain characteristics. There are 3 methods available: Inflation-normalization ('IN'), Bistochastic Scaling ('BS'), and Stationary Distribution Scaling ('SDS').
#'
#' See the [manuscript](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10324253/) for more details on the AMEND algorithm. There have been modifications to the algorithm as presented in the original paper.
#'
#' @param graph,adj_matrix,edge_list A single graph like object in the form of an igraph, adjacency matrix, or edge list. Or, a named list containing multiple graph-like objects of the same type, to be merged. The merged graph must be connected. If a third column is provided in an edge list, these are taken as edge weights. Only one of graph, adj_matrix, or edge_list should be given, with priority given to graph, then adj_matrix, then edge_list. See 'Details' for correct naming conventions.
#' @param n Size of the final module to be approximated.
#' @param data A named list of named numeric vectors (list elements correspond to graph components), a named numeric vector, or a character scalar (denoting a vertex attribute of the input igraph object) containing the experimental data from which seed values for RWR will be derived according to _FUN_ and _FUN.params_ args. See 'Details' for correct naming conventions.
#' @param node_type A named list of named character vectors (list elements correspond to graph components), a named character vector, a character scalar (denoting a vertex attribute of the input igraph object), or NULL. Denotes the component and the multiplex layer of each node. If NULL and multiplex and/or heterogeneous, node labels must follow 'name|type_layer' naming scheme (e.g., MYC|gene_1). See 'Details' for correct naming conventions.
#' @param brw.attr A named list of named numeric vectors (list elements correspond to graph components), a named numeric vector, a character vector (of nodes of interest), a character scalar (denoting a vertex attribute of the input igraph object), or NULL. Biased random walk vertex attribute values should be non-negative, with larger values increasing the transition probabilities to a node in RWR. If NULL, all nodes are given a value of 1. Character vector gives nodes of interest towards which transition probabilities are increased. See 'Details' for biased random walk info.
#' @param FUN A function, named list of functions, named list of character scalars, a single character scalar, or NULL. Function for transforming values in _data_ to derive seed values for RWR. Names correspond to graph components. Character strings should correspond to a default function: one of 'binary', 'shift_scale', 'p_value', or 'exp'. NULL means no transformation is done to values in _data_. See 'Details' for descriptions of default functions.
#' @param FUN.params A named list of lists of named function arguments, a named list of named function arguments, or NULL. Function arguments to be passed to _FUN_. Names should match names in _FUN_.
#' @param heterogeneous Logical. If TRUE, graph is considered heterogeneous (more than one distinct node type, e.g., proteins and metabolites), and _node_type_ must be included as an argument or graph vertex attribute.
#' @param multiplex Logical. If true, graph is assumed to contain multiplex components.
#' @param aggregate.multiplex A named list (for aggregating) or NULL (for no aggregating). The list element 'primary' contains the name of the primary layer for a multiplex component whose edges will be used during subnetwork identification. The multiplex component is collapsed onto this primary layer. The list element 'agg.method' contains a character scalar referring to an aggregation function for aggregating vertex attributes ('mean', 'median', 'sum', 'gmean', 'hmean'). gmean=geometric, hmean=harmonic.
#' @param normalize Normalization scheme of adjacency matrix for random walk with restart
#' @param k Value between 0 and 1. When normalize = "modified_degree", the adjacency matrix is first left and right multiplied by a diagonal matrix of node degrees, which is raised to the power -_k_. As _k_ increases, edge weights are penalized more for the degrees of their adjacent nodes.
#' @param degree.bias List or NULL. NULL (default) for no degree bias adjustment. List names should include 'component', which is a character string of length greater than or equal to 1 specifying the components or layers to which the degree bias adjustment will be applied, and 'method', which is a character string of 'BS' (bistochastic scaling), 'IN' (inflation-normalization), or 'SDS' (stationary distribution scaling).
#' @param jump.prob,net.weight A named vector, or NULL. _jump.prob_ is the probability of random walker jumping from one component of graph to another in RWR. _net.weight_ is the relative weight given to nodes of a component of graph, applied to seed vector in RWR. Only used when heterogeneous=TRUE.
#' @param switch.layer.prob,layer.weight A named list of named vectors, or NULL. _switch.layer.prob_ is the probability of random walker to switch from current layer in a multiplex to another layer in same component. _layer.weight_ is relative weight given to nodes of a layer of a component of graph, applied to seed vector in RWR. List element names correspond to multiplex components, and vector names correspond to layers within a multiplex.
#' @param verbose Logical. Whether to output current iteration number to show progress.
#' @param eta Starting filtering rate. If NULL (default), a value is chosen based on the input network size and parameter _n_.
#' @param identifier Optional. For use when performing many runs of AMEND to keep track of progress.
#' @param in.parallel Logical. Run computations in parallel.
#' @param n.cores number of cores to use in parallel computations
#'
#' @return a named list with the following elements:
#' * module: the final module (i.e., subnetwork)
#' * score: final module score
#' * subnetworks: a list of node names contained in intermediate subnetworks
#' * stats: network statistics
#' * time: run time
#' * input_params: list of input parameters
#'
#' @seealso [create_integrated_graph()], [list_input_checks()], [transition_matrix()], [RWR()], [heinz()], [get_subnetwork()], [inflate_normalize()], [bistochastic_scaling()]
#'
#' @examples
#' # Attach igraph library
#' library(igraph)
#'
#' # Inspect the igraph object included in AMEND package.
#' # One can see it has vertex attributes name, symbol, and ECI
#' glut4_graph
#' head(V(glut4_graph)$ECI)
#'
#' \dontrun{
#' # Use run_AMEND() with an igraph object with a vertex attribute matching data arg
#' subnet1 = run_AMEND(graph = glut4_graph, n = 25, data = "ECI", FUN = "shift_scale",
#'                     FUN.params = list(DOI = -1, w = 0.5), normalize = "degree")
#'
#' # Use run_AMEND() with an adjacency matrix and a vector of node scores
#' subnet2 = run_AMEND(adj_matrix = glut4_adjM, data = eci_scores, FUN = "shift_scale",
#'                     FUN.params = list(DOI = -1, w = 0.5), normalize = "degree")
#'
#' #=== Monoplex-heterogeneous network ===#
#' # Using a list of edge lists, with data as list input
#' # Create object for 'data' argument
#' data.list.input = vector("list", 2)
#' names(data.list.input) = c("prot", "meta")
#' for(i in seq_along(data.list.input)){
#'   uniq.names = extract_string(unique(c(edgelists[[i]][,1], edgelists[[i]][,2])), "\\|", 1)
#'   data.list.input[[i]] = rnorm(n = length(uniq.names), mean = 0, sd = 1)
#'   names(data.list.input[[i]]) = uniq.names
#' }
#'
#' # Setting multiplex/heterogeneous parameters
#' net.weight = c(prot = 0.5, meta = 0.5)
#' jump.prob = c(prot = 0.7, meta = 0.3)
#'
#' # FUN and FUN.params
#' FUN = list(prot = "shift_scale", meta = function(x, k) abs(x) * k)
#' FUN.params = list(prot = list(DOI = 1, w = 0.5), meta = list(k = 1.2))
#' subnet3 = run_AMEND(edge_list = edgelists, n = 50, data = data.list.input, node_type = NULL,
#'                     brw.attr = NULL, FUN = FUN, FUN.params = FUN.params, heterogeneous = TRUE,
#'                     multiplex = FALSE, normalize = "degree", jump.prob = jump.prob,
#'                     net.weight = net.weight, verbose = TRUE)
#'
#' #=== Multiplex-heterogeneous network ===#
#' # List of graphs, data as vertex attribute
#' # Create data vertex attributes, called 'scores'
#' v.attr.name = "scores"
#' node.types = c("prot_1", "prot_2", "prot_3", "meta")
#' for(i in seq_along(node.types)){
#'   igraph::vertex_attr(list.of.graphs[[node.types[i]]],
#'                       v.attr.name) = runif(vcount(list.of.graphs[[node.types[i]]]))
#' }
#'
#' # Setting multiplex/heterogeneous parameters
#' layer.weight = list(prot = rep(1/3, 3))
#' switch.layer.prob = list(prot = c(prot_1 = 0.5, prot_2 = 0.2, prot_3 = 0.8))
#' net.weight = c(prot = 0.5, meta = 0.5)
#' jump.prob = c(prot = 0.5, meta = 0.3)
#'
#' subnet4 = run_AMEND(graph = list.of.graphs, n = 100, data = "scores", node_type = NULL,
#'                     brw.attr = NULL, FUN = "p_value", FUN.params = NULL,
#'                     heterogeneous = TRUE, multiplex = TRUE, normalize = "modified_degree",
#'                     k = 0.5, jump.prob = jump.prob, net.weight = net.weight,
#'                     switch.layer.prob = switch.layer.prob,
#'                     layer.weight = layer.weight, verbose = TRUE)
#'
#' #=== Multiplex-homogeneous network ===#
#' # List of adjacency matrices, data as list input
#' # Create a list of binary vectors (e.g., significantly up- or down-regulated genes)
#' data.list.input = vector("list", 3)
#' names(data.list.input) = c("prot_1", "prot_2", "prot_3")
#' for(i in seq_along(data.list.input)){
#'   data.list.input[[i]] = sample(x = c(0,1), size = nrow(list.of.adjmats[[i]]),
#'                                 replace = T, prob = c(0.9, 0.1)) # binary
#'   names(data.list.input[[i]]) = rownames(list.of.adjmats[[i]])
#' }
#'
#' # Setting multiplex/heterogeneous parameters
#' layer.weight = list(prot = c(prot_1 = 0.4, prot_2 = 0.3, prot_3 = 0.3))
#' switch.layer.prob = list(prot = c(prot_1 = 0.4, prot_2 = 0.5, prot_3 = 0.5))
#'
#' subnet5 = run_AMEND(adj_matrix = list.of.adjmats, n = 100, data = data.list.input,
#'                     node_type = NULL, brw.attr = NULL, FUN = "binary",
#'                     FUN.params = NULL, heterogeneous = FALSE, multiplex = TRUE,
#'                     normalize = "modified_degree", k = 0.5, layer.weight = layer.weight,
#'                     switch.layer.prob = switch.layer.prob, verbose = TRUE)
#' }
#'
#' @export
run_AMEND <- function(graph = NULL, adj_matrix = NULL, edge_list = NULL, n = 50, data = NULL, node_type = NULL, brw.attr = NULL,
                      FUN = NULL, FUN.params = NULL, heterogeneous = FALSE, multiplex = FALSE, aggregate.multiplex = NULL,
                      normalize = c("degree", "modified_degree"), k = 0.5, degree.bias = NULL,
                      jump.prob = NULL, net.weight = NULL, switch.layer.prob = NULL, layer.weight = NULL,
                      verbose = FALSE, eta = NULL, identifier = 1, in.parallel = FALSE, n.cores){
  start_time = Sys.time()

  # Variables not set by user, but could change in future
  ma_window = 5 # Moving average window... set to Inf for a running average

  normalize = match.arg(normalize)

  # Handling brw.attr argument
  if(is.character(brw.attr) && length(brw.attr) == 1){
    brw.attr.nm = brw.attr
  }else brw.attr.nm = "brw.values"

  # Verifying that degree.bias argument is a list with at least a 'method' named element
  if(!is.null(degree.bias) && (!is.list(degree.bias) || !'method' %in% names(degree.bias))) stop('\'degree.bias\' should be a list with names \'method\' and \'component\'.\n\'method\' should be one of \'BS\', \'IN\', or \'SDS\'.')

  # Creating integrated graph
  graph = create_integrated_graph(graph = graph, adj_matrix = adj_matrix, edge_list = edge_list, data = data, node_type = node_type, brw.attr = brw.attr,
                                  FUN = FUN, FUN.params = FUN.params, heterogeneous = heterogeneous, multiplex = multiplex, lcc = TRUE)

  if(is.character(brw.attr)){
    name.only = get.type(V(graph)$name, 1)
    name.comp = paste(name.only, get.type(V(graph)$name, 3), sep = "|")
    name.comp.layer = V(graph)$name
    if(any(brw.attr %in% c(name.only, name.comp, name.comp.layer))) brw.flag = TRUE else brw.flag = FALSE
  }else brw.flag = FALSE
  graph = set_brw_attr(graph = graph, brw.attr = brw.attr, brw.attr.nm = brw.attr.nm, brw.flag = brw.flag)

  # Sort order of vertices in graph
  nt = unique(sort(V(graph)$node_type))
  tmp = c()
  for(i in seq_along(nt)) tmp = c(tmp, which(V(graph)$node_type == nt[i]))
  p = order(tmp)
  graph = igraph::permute(graph, p)
  comps.orig = unique(get.type(V(graph)$name, 3))

  # Create aggregated graph
  if(is.null(aggregate.multiplex)){
    agg.graph = graph
  }else{
    # Check for correct entries in aggregate.multiplex
    if(is.character(aggregate.multiplex)){
      aggregate.multiplx = list(primary = aggregate.multiplex, agg.method = "mean")
    }
    if(is.list(aggregate.multiplex) && !is.null(names(aggregate.multiplex))){
      if(!"primary" %in% names(aggregate.multiplex)) stop("Names of aggregate.multiplex must include 'primary' and 'agg.method'.")
      if(!"agg.method" %in% names(aggregate.multiplex)) aggregate.multiplex$agg.method = "mean"
      if(length(aggregate.multiplex$agg.method) > 1) aggregate.multiplex$agg.method = aggregate.multiplex$agg.method[1]
      if(!any(aggregate.multiplex$agg.method %in% c("mean", "median", "sum", "gmean", "hmean"))){
        aggregate.multiplex$agg.method = "mean"
        warning("Unrecognized aggregation function name. Using mean.\nMust be one of c('mean', 'median', 'sum', 'gmean', 'hmean').")
      }
      if(!any(aggregate.multiplex$primary %in% V(graph)$node_type)) stop("Unrecognized primary 'component_layer' for multiplex aggregation.")
      # Making sure only one layer per component is given as primary
      aggregate.multiplex$primary = unique(aggregate.multiplex$primary)
      tbl.tmp = table(extract_string(aggregate.multiplex$primary, "_", 1))
      if(any(tbl.tmp > 1)){
        warning("More than one primary layer given for multiplex aggregation. Choosing first layer given.")
        id = match(names(tbl.tmp), extract_string(aggregate.multiplex$primary, "_", 1))
        aggregate.multiplex$primary = aggregate.multiplex$primary[id]
      }
    }else stop("Incorrect input type for aggregate.multiplex.")
    agg.graph = create_aggregated_graph(graph = graph, control = aggregate.multiplex)
    # Create new vertex attr that contains name|component info of aggregated multiplex components
    for(i in seq_along(aggregate.multiplex$primary)){
      comp = extract_string(aggregate.multiplex$primary[i], "_", 1)
      tmp.nm = paste(get.type(V(graph)$name, 1), comp, sep = "|")
      tmp.id = which(get.type(V(graph)$name, 3) == comp)
      igraph::vertex_attr(graph, "agg_name", tmp.id) = tmp.nm[tmp.id]
      igraph::vertex_attr(graph, "agg_name") = ifelse(is.na(igraph::vertex_attr(graph, "agg_name")), V(graph)$name, igraph::vertex_attr(graph, "agg_name"))
    }
  }

  # Get adjacency matrix
  adj_mat = igraph::as_adjacency_matrix(graph = graph, attr = "weight", sparse = TRUE)

  # Additional checks for multiplex and heterogeneous parameters.
  if(multiplex){
    graph_layer = unique(V(graph)$node_type)
    graph_type = extract_string(graph_layer, "_", 1) # node types (NOT layer specific)
    multiplex.comps = unique(graph_type[grepl("_", graph_layer)])

    ## Checks for 'layer.weight' object (analogous to jump.prob)
    # Values for each multiplex component must be non-negative and sum to one
    # layer.weight, multiplex.comps, names(graph) (where graph is a list)
    if(is.null(layer.weight)){
      layer.weight = vector("list", length(multiplex.comps))
      names(layer.weight) = multiplex.comps
      for(i in seq_along(layer.weight)){
        nL = sum(graph_type == multiplex.comps[i])
        layer.weight[[i]] = rep(1/nL, nL)
        names(layer.weight[[i]]) = graph_layer[graph_type == multiplex.comps[i]]
      }
    }else if(is.list(layer.weight) && !is.null(names(layer.weight))){
      for(i in seq_along(layer.weight)){
        if(!names(layer.weight)[i] %in% multiplex.comps) stop(paste0("'", names(layer.weight)[i], "' from 'layer.weight' does not match any of the multiplex graph types given: ", paste(multiplex.comps, collapse = ", ")))
        nL = sum(graph_type == names(layer.weight)[i])
        if(length(layer.weight[[i]]) != nL){
          stop("Incorrect number of elements in layer.weight[[", i, "]]. Should be ", nL, ".")
        }else if(is.null(names(layer.weight[[i]]))){
          warning("No layer names given in layer.weight[[", i, "]]. Assuming layer order matches order of input for graph/adj_matrix/edge_list.")
          names(layer.weight[[i]]) = graph_layer[graph_type == names(layer.weight)[i]]
        }else if(any(!names(layer.weight[[i]]) %in% c(graph_layer, extract_string(graph_layer[grep("_", graph_layer)], "_", 2)))) stop("Names in layer.weight[[", i, "]] don't match any names in graph.")
        if(any(!names(layer.weight[[i]]) %in% graph_layer & names(layer.weight[[i]]) %in% extract_string(graph_layer[grep("_", graph_layer)], "_", 2))){
          id = which(!names(layer.weight[[i]]) %in% graph_layer & names(layer.weight[[i]]) %in% extract_string(graph_layer[grep("_", graph_layer)], "_", 2))
          names(layer.weight[[i]])[id] = paste(names(layer.weight)[i], names(layer.weight[[i]])[id], sep='_')
        }
        # Check values
        if(sum(layer.weight[[i]]) != 1) stop("Elements of layer.weight[[", i, "]] do not sum to 1.")
        if(any(layer.weight[[i]] < 0)) stop("Elements of layer.weight[[", i, "]] are not non-negative.")
      }
      if(any(!multiplex.comps %in% names(layer.weight))){ # If any multiplex components are missing
        id = which(!multiplex.comps %in% names(layer.weight))
        warning(paste0("No 'layer.weight' values given for multiplex graph type(s): ", paste(multiplex.comps[id], collapse = ", "), ". Defaulting to uniform weights."))
        for(i in id){
          nL = sum(graph_type == multiplex.comps[i])
          layer.weight[[length(layer.weight) + 1]] = rep(1/nL, nL)
          names(layer.weight[[length(layer.weight)]]) = graph_layer[graph_type == multiplex.comps[i]]
          names(layer.weight)[length(layer.weight)] = multiplex.comps[i]
        }
      }
    }else stop("layer.weight must be a named list or NULL.")

    ## Checks for 'switch.layer.prob' object (analogous to jump.prob). Values must be between 0 and 1.
    # Need a value for each layer of each multiplex component. They do not need to sum to one within a multiplex. Represent probability of switching from current layer to any other layer.
    # NULL or named list
    # If NULL, assign 0.5 to all layers in each multiplex component
    # If named list, check that all names match multiplex components and check that all layer names match.
    # switch.layer.prob, multiplex.comps, names(graph) (where graph is a list)
    if(is.null(switch.layer.prob)){
      # Assign 0.5 to all layers in each multiplex component
      switch.layer.prob = vector("list", length(multiplex.comps))
      names(switch.layer.prob) = multiplex.comps
      for(i in seq_along(switch.layer.prob)){
        nL = sum(graph_type == multiplex.comps[i])
        switch.layer.prob[[i]] = rep(0.5, nL)
        names(switch.layer.prob[[i]]) = graph_layer[graph_type == multiplex.comps[i]]
      }
    }else if(is.list(switch.layer.prob) && !is.null(names(switch.layer.prob))){
      for(i in seq_along(switch.layer.prob)){
        if(!names(switch.layer.prob)[i] %in% multiplex.comps) stop(paste0("Names of switch.layer.prob do not match any of the multiplex graph types given: ", paste(multiplex.comps, collapse = ", ")))
        nL = sum(graph_type == names(switch.layer.prob)[i])
        if(length(switch.layer.prob[[i]]) == 1){
          # Apply this value to all layers
          switch.layer.prob[[i]][2:nL] = switch.layer.prob[[i]][1]
          names(switch.layer.prob[[i]]) = graph_layer[graph_type == names(switch.layer.prob)[i]]
        }else if(is.null(names(switch.layer.prob[[i]]))){
          if(length(switch.layer.prob[[i]]) != nL) stop("Incorrect number of layers given in switch.layer.prob[[", i, "]]")
          warning("No layer names given in switch.layer.prob[[", i, "]]. Assuming layer order matches order of input for graph/adj_matrix/edge_list.")
          names(switch.layer.prob[[i]]) = graph_layer[graph_type == names(switch.layer.prob)[i]]
        }else if(any(!names(switch.layer.prob[[i]]) %in% c(graph_layer, extract_string(graph_layer[grep("_", graph_layer)], "_", 2)))) stop("Names in switch.layer.prob[[", i, "]] don't match any names in graph.")
        if(any(!names(switch.layer.prob[[i]]) %in% graph_layer & names(switch.layer.prob[[i]]) %in% extract_string(graph_layer[grep("_", graph_layer)], "_", 2))){
          id = which(!names(switch.layer.prob[[i]]) %in% graph_layer & names(switch.layer.prob[[i]]) %in% extract_string(graph_layer[grep("_", graph_layer)], "_", 2))
          names(switch.layer.prob[[i]])[id] = paste(names(switch.layer.prob)[i], names(switch.layer.prob[[i]])[id], sep='_')
        }
        # Check range
        if(any(switch.layer.prob[[i]] < 0) || any(switch.layer.prob[[i]] > 1)) stop("Incorrect range for values of switch.layer.prob[[", i, "]]. Must be between 0 and 1.")
      }
      if(any(!multiplex.comps %in% names(switch.layer.prob))){ # If any multiplex components are missing
        id = which(!multiplex.comps %in% names(switch.layer.prob))
        warning(paste0("No 'switch.layer.prob' values given for multiplex graph type(s): ", paste(multiplex.comps[id], collapse = ", "), ". Defaulting to 0.5."))
        for(i in id){
          nL = sum(graph_type == multiplex.comps[i])
          switch.layer.prob[[length(switch.layer.prob) + 1]] = rep(0.5, nL)
          names(switch.layer.prob[[length(switch.layer.prob)]]) = graph_layer[graph_type == multiplex.comps[i]]
          names(switch.layer.prob)[length(switch.layer.prob)] = multiplex.comps[i]
        }
      }
    }else stop("switch.layer.prob must be a named list or NULL.")
  }
  if(heterogeneous){
    uniq.types = unique(extract_string(unique(V(graph)$node_type), "_", 1))

    # Checks for jump.prob, used in transition_matrix() and RWR()
    # Must be a named numeric vector or NULL.
    # Each element in jump.prob should be a number between 0 and 1
    # If NULL, it defaults to 0.5 for each node type
    # If vector, name needs to match to name in names(graph)
    if(is.null(jump.prob)){
      jump.prob = rep(0.5, length(uniq.types))
      names(jump.prob) = uniq.types
    }else if(is.numeric(jump.prob) && !is.null(names(jump.prob))){
      if(any(!names(jump.prob) %in% uniq.types)) stop("Names of 'jump.prob' don't match names of graph.")
      # Check that values are between 0 and 1
      for(i in seq_along(jump.prob)) if(any(jump.prob[i] > 1) || any(jump.prob[i] < 0)) stop("Incorrect range for elements in jump.prob[", i, "]. Must be between 0 and 1.")
      if(any(!uniq.types %in% names(jump.prob))){ # If any node types are missing
        id = which(!uniq.types %in% names(jump.prob))
        warning(paste0("No 'jump.prob' values given for node type(s) ", paste(uniq.types[id], collapse = ", "), ". Defaulting to 0.5."))
        for(i in id){
          jump.prob[length(jump.prob) + 1] = 0.5
          names(jump.prob)[length(jump.prob)] = uniq.types[i]
        }
      }
    }else stop("Incorrect input type for 'jump.prob'. Must be a named numeric vector or NULL.")

    # Check for 'net.weight', weights to give seed values of each node type for RWR
    # named numeric vector or NULL
    # If NULL, give each node type (not considering layers of multiplexes) equal weight
    # Values must sum to one and be non-negative
    # Names of vector must match names of graph
    if(is.null(net.weight)){
      net.weight = rep(1 / length(uniq.types), length(uniq.types))
      names(net.weight) = uniq.types
    }else if(is.numeric(net.weight) && !is.null(names(net.weight))){
      if(!all(names(net.weight) %in% uniq.types)) stop("'net.weight' has node types not contained in graph.")
      if(!all(uniq.types %in% names(net.weight))) stop("Not all node types are in 'net.weight'.")
      if(sum(net.weight) != 1) stop("Elements of 'net.weight' do not sum to 1.")
      if(any(net.weight < 0)) stop("Elements of 'net.weight' are not non-negative.")
    }else stop("Incorrect input type for 'net.weight'. Must be a named numeric vector or NULL.")
  }
  repeat{
    all_scores <- list()
    all_nets <- list()

    subg <- vector(mode = "list", length = 2)
    subg[[1]] <- graph
    sub.ag = vector("list", length = 2)
    sub.ag[[1]] = agg.graph

    iter.num <- 1

    avg.rd = 0.07 # initial average rate difference

    if(is.null(eta)) eta = get.eta0(agg.graph, n, "log.ratio")

    if(verbose) message(paste("Starting filtering rate:", round(eta, 4)))
    repeat{
      # if(iter.num == 1) browser()
      if(verbose) message(paste("Iteration:", iter.num))
      if(iter.num == 1){
        rate.difference = avg.rd
        e = eta
        decay.N = vcount(agg.graph)
        adj <- adj_mat
      }else{
        rate.difference = c(avg.rd, unlist(lapply(all_scores, function(x) x[10] - x[9])))
        rate.difference = rate.difference[ifelse(length(rate.difference) >= ma_window, length(rate.difference) - ma_window + 1, 1):length(rate.difference)]
        e = all_scores[[iter.num-1]][9]

        if(iter.num == 2){
          decay.N = vcount(agg.graph)
        }else decay.N = length(all_nets[[iter.num-2]])

        adj <- adj[-rm.id, -rm.id, drop=FALSE]
      }
      # Setting shifting percentile (i.e., filtering rate) for this iteration
      decay <- find_decay(N = decay.N, eta0 = e, n = n, rate.diff = rate.difference)
      l <- exp_filtering_rate(eta0 = e, d = decay)
      ll <- l[ifelse(iter.num == 1, 1, 2)]

      # !!!! CHANGE THESE 'IF' STATEMENTS TO CHECK IF A COMPONENT OR LAYER HAS BEEN COMPLETELY REMOVED !!!!
      # Determining whether graph is still heterogeneous, b/c all of one node type may have been filtered out
      if(heterogeneous){
        comps = unique(get.type(V(subg[[1]])$name, 3))
        if(length(comps) < 2){ # length(unique(V(subg[[1]])$node_type)) < 2
          heterogeneous = FALSE
          if(verbose) message(paste0("All nodes of type \'", setdiff(comps.orig, comps), "\' have been filtered out."))
        }
      }
      # Determine whether graph is still multiplex. Check the number of layers for each node type
      if(multiplex){
        comps = unique(get.type(V(subg[[1]])$name, 3))
        tmp = FALSE
        for(i in seq_along(comps)){
          id = which(get.type(V(subg[[1]])$name, 3) == comps[i])
          if(length(unique(V(subg[[1]])$node_type[id])) > 1) tmp = TRUE
        }
        multiplex = tmp
        if(!multiplex) message("Graph is now monoplex (i.e., all node types have only one layer).")
      }

      # Normalize adjacency matrix to get transition matrix
      n.adjM = transition_matrix(adjM = adj, norm = normalize, k = k, heterogeneous = heterogeneous, multiplex = multiplex,
                                 brw.attr = igraph::vertex_attr(subg[[1]], brw.attr.nm),
                                 jump.prob = jump.prob, switch.layer.prob = switch.layer.prob,
                                 degree.bias = degree.bias, in.parallel = in.parallel, n.cores = n.cores)

      # Create matrix of seed values
      Seeds = matrix(V(subg[[1]])$seeds, ncol = 1, dimnames = list(V(subg[[1]])$name))

      # Choosing Restart parameter value through a 'normal' grid search
      rgs = restart_grid_search(ig = subg[[1]], ag = sub.ag[[1]], n.adj.M = n.adjM, seeds = Seeds, filtering_rate = ll, heterogeneous = heterogeneous,
                                multiplex = multiplex, net.weight = net.weight, layer.weight = layer.weight, iteration = iter.num, agg.method = aggregate.multiplex, degree.bias = degree.bias, in.parallel = in.parallel, n.cores = n.cores)
      subg[[2]] = expand_graph(ig = subg[[1]], ag = rgs[[1]], control = aggregate.multiplex)
      sub.ag[[2]] = rgs[[1]]

      # Update brw attr for full graph
      subg[[2]] = set_brw_attr(graph = subg[[2]], brw.attr = brw.attr, brw.attr.nm = brw.attr.nm, brw.flag = brw.flag)

      # break out of repeat loop if there is no change in network or if subnet size <= 2
      # if(vcount(subg[[2]]) == vcount(subg[[1]]) || vcount(subg[[2]]) <= 2) break
      if(vcount(sub.ag[[2]]) == vcount(sub.ag[[1]]) || vcount(sub.ag[[2]]) <= 2) break

      mCC = rgs[[2]][1] # Mean core-clustering coefficient of aggregated graph
      mZ = rgs[[2]][2] # Mean standardized experimental scores of aggregated graph
      n_nodes = rgs[[2]][3] # Size of subnetwork from aggregated graph

      sn_score = mCC * mZ # subnetwork score for aggregated graph

      all_scores[[length(all_scores) + 1]] = c(decay, rgs[[3]], sn_score, mZ, mCC, n_nodes, igraph::ecount(sub.ag[[2]]), igraph::edge_density(sub.ag[[2]]), ll, mean(!V(sub.ag[[1]])$name %in% V(sub.ag[[2]])$name))
      # all_nets[[length(all_nets) + 1]] = V(subg[[2]])$name
      all_nets[[length(all_nets) + 1]] = V(sub.ag[[2]])$name

      # IDs of nodes that were removed this iteration
      rm.id = which(!V(subg[[1]])$name %in% V(subg[[2]])$name)

      iter.num = iter.num + 1
      subg[[1]] = subg[[2]]
      sub.ag[[1]] = sub.ag[[2]]
    }
    all_scores <- do.call("rbind", all_scores)
    colnames(all_scores) <- c("Decay", "Restart parameter", "Network score", "Avg Z", "Avg CCC", "Nodes", "Edges", "Density",
                              "Filtering rate", "Observed filtering rate")
    all_scores <- round(as.data.frame(all_scores), 3)

    # final module cannot be larger than twice the user-specified 'n'
    size.cond = all_scores[,6] <= 2 * n
    if(any(size.cond)){
      break
    }else{
      if(verbose) message(paste0("Algorithm didn't converge with eta=", round(eta, 4), ". Trying a larger starting filtering rate. ID=", identifier))
      if(eta == 1) stop("Algorithm is unable to converge. Try a larger 'n'.")
      eta = eta * 1.1 # increase by 10%
      if(eta > 1) eta = 1
    }
  }

  # Choosing subnetwork with maximum score and (in the case of ties) closest to final module size 'n'
  # and in case of two networks with same max score and same distance from 'n', choose larger subnetwork
  if(sum(all_scores[size.cond, 3] == max(all_scores[size.cond, 3])) > 1){ # Handling ties
    d = abs(all_scores[,6] - n)
    ids = which(all_scores[,3] == max(all_scores[size.cond, 3]) & size.cond)
    max.id = which(all_scores[,3] == max(all_scores[size.cond, 3]) & d == min(d[ids]) & size.cond)[1]
  }else max.id = which(all_scores[,3] == max(all_scores[size.cond, 3]) & size.cond)

  # best_sn = igraph::induced_subgraph(graph, which(V(graph)$name %in% all_nets[[max.id]]))
  best_sn = igraph::induced_subgraph(agg.graph, which(V(agg.graph)$name %in% all_nets[[max.id]]))
  best_score = all_scores[max.id,3]

  # Remove text after '|' for simple graphs
  if(!multiplex && !heterogeneous){
    V(best_sn)$name = extract_string(V(best_sn)$name, "\\|", 1)
    all_nets = lapply(all_nets, function(y) extract_string(y, "\\|", 1))
  }

  end_time <- Sys.time()
  time <- end_time - start_time
  if(verbose){
    message(paste0("*** Converged! *** ID=", identifier))
  }else message("*** Converged! ***")

  input_params = list(n = n, normalize = normalize, k = k, FUN = FUN, FUN.params = FUN.params, heterogeneous = heterogeneous, multiplex = multiplex, jump.prob = jump.prob, net.weight = net.weight,
                      switch.layer.prob = switch.layer.prob, layer.weight = layer.weight)

  return(list(module = best_sn, score = best_score, subnetworks = all_nets, stats = all_scores, time = time, input_params = input_params))
}

#' @title Grid search for the restart parameter in RWR
#'
#' @description
#' Grid search for the restart parameter in RWR. For use inside `run_AMEND()`
#'
#' @details
#' This is a simple grid search for the restart probability parameter in RWR. For a grid of restart probability vaules, `RWR()` is run, raw scores are shifted by the current filtering rate of that iteration, then `heinz()` is run to find optimal subnetworks, which are scored. The restart probability resulting in largest subnetwork score is used.
#'
#' @param ig Input graph
#' @param ag Aggregated graph
#' @param n.adj.M Normalized adjacency matrix
#' @param seeds A named vector of seed values
#' @param filtering_rate Quantile for shifting the raw RWR scores
#' @param heterogeneous Logical. If TRUE, network is considered heterogeneous (two distinct node types, e.g., proteins and metabolites). If TRUE, _node_type_ must be included as an argument or graph vertex attribute
#' @param multiplex Logical. If true, graph is assumed to contain multiplex components.
#' @param net.weight A named vector, or NULL. Relative weight given to nodes of a component of graph, applied to seed vector in RWR. Only used when heterogeneous=TRUE.
#' @param layer.weight A named list of named vectors, or NULL. Relative weight given to nodes of a layer of a component of graph, applied to seed vector in RWR. List element names correspond to multiplex components, and vector names correspond to layers within a multiplex.
#' @param iteration Current iteration of the AMEND algorithm.
#' @param agg.method A named list. The element 'primary' contains the name of the primary layer for a multiplex component to be used during subnetwork identification. The element 'agg.method' contains a character scalar referring to an aggregation function.
#' @param degree.bias List or NULL. NULL (default) for no degree bias adjustment. List names should include 'component', which is a character string of length greater than or equal to 1 specifying the components or layers to which the degree bias adjustment will be applied, and 'method', which is a character string of 'BS' (bistochastic scaling), 'IN' (inflation-normalization), or 'SDS' (stationary distribution scaling).
#' @param in.parallel Logical. Run computations in parallel.
#' @param n.cores number of cores to use in parallel computations
#'
#' @return a list containing an igraph object, subnetwork score, and restart value
#'
#' @examples
#' library(igraph)
#'
#' # Creating an adjacency matrix and graph
#' adjm = matrix(c(0, 1, 0, 0, 0, 0, 0, 0,
#'                1, 0, 1, 1, 1, 0, 0, 0,
#'                0, 1, 0, 0, 1, 1, 0, 0,
#'                0, 1, 0, 0, 0, 0, 1, 0,
#'                0, 1, 1, 0, 0, 0, 1, 1,
#'                0, 0, 1, 0, 0, 0, 0, 0,
#'                0, 0, 0, 1, 1, 0, 0, 0,
#'                0, 0, 0, 0, 1, 0, 0, 0), nrow = 8, dimnames = list(1:8, 1:8))
#' g = graph_from_adjacency_matrix(adjm, mode = 'undirected')
#' V(g)$name = 1:8
#' plot(g)
#'
#' # Normalizing the adjacency matrix
#' adj_norm = adjm %*% diag(1 / degree(g))
#' dimnames(adj_norm) = dimnames(adjm)
#'
#' # Creating a named vector of seed values
#' seeds = runif(8)
#' names(seeds) = 1:8
#'
#' search = AMEND:::restart_grid_search(ig = g, ag = g, n.adj.M = adj_norm,
#'                                      seeds = seeds, filtering_rate = 0.3)
#' plot(search[[1]])
#'
restart_grid_search <- function(ig, ag, n.adj.M, seeds, filtering_rate, heterogeneous = FALSE, multiplex = FALSE, net.weight, layer.weight, iteration = 1, agg.method = NULL, degree.bias = NULL, in.parallel = FALSE, n.cores){
  if(!is.null(agg.method)){
    agg.fun = get_aggregate_method(agg.method$agg.method)
    flag = TRUE
  }else flag = FALSE

  if(iteration == 1){
    grid <- seq(0.5, 0.95, by = 0.05)
  }else grid <- seq(0.2, 0.95, by = 0.05)

  # Degree bias adjustment
  if(!is.null(degree.bias) && degree.bias$method == 'SDS'){
    sds = RWR(nadjM = n.adj.M, setSeeds = seeds, restart = 0, heterogeneous = heterogeneous, multiplex = multiplex, net.weight = net.weight, layer.weight = layer.weight)
    sds = as.vector(sds)
    names(sds) = V(ig)$name
    sds[sds == 0] = min(sds[sds != 0]) # To avoid dividing by zero
  }

  if(0){
    #=======================================================================#
    #=== Stouffer's Z for aggregating centrality measures (Experimental) ===#
    # 0. Get input arguments
    gamma = 0.75
    w = c(ccc = 1, hc = 1, bc = 1)
    # 1. Calculate core-clustering coefficient, harmonic centrality, and betweenness for each node
    # cutoff could be set to a certain proportion of the diameter of the graph (i.e., maximum eccentricity, i.e., maximum of the maximum shortest path length for all nodes)
    g.d = igraph::diameter(graph = ag)
    kpath = ifelse(gamma == 1, -1, gamma * g.d)
    ccc = core_cc(ag)
    hc = igraph::harmonic_centrality(graph = ag, mode = 'out', normalized = FALSE, cutoff = kpath)
    bc = igraph::betweenness(graph = ag, directed = TRUE, normalized = FALSE, cutoff = kpath)
    centrality.matrix = matrix(c(ccc, hc, bc), ncol = 3)
    # 2. Get empirical cumulative probabilities from these values.
    centrality.matrix = apply(centrality.matrix, 2, function(x) stats::ecdf(x)(x))
    # 3. Convert to standard normal values with inverse standard normal CDF
    centrality.matrix = apply(centrality.matrix, 2, stats::qnorm)
    # 4. For each node, compute Stouffer's Z-score by taking a weighted mean of these standard normal values.
    aggregated.centrality = apply(centrality.matrix, 1, function(x) sum(x * w[match(names(x), names(w))]) / sqrt(sum(w^2)))
    ## Note: Allow the user to supply positive weights for each centrality measure.
    #=============================#
  }

  if(!in.parallel){
    score_metrics <- matrix(nrow = length(grid), ncol = 4)
    nets <- vector(mode = "list", length = length(grid))
    rm.id = NULL
    for(i in seq_along(grid)){
      # RWR
      PTmatrix2 <- RWR(nadjM = n.adj.M, setSeeds = seeds, restart = grid[i], heterogeneous = heterogeneous, multiplex = multiplex, net.weight = net.weight, layer.weight = layer.weight)
      rwr_score <- as.vector(PTmatrix2)
      names(rwr_score) <- V(ig)$name

      # Degree bias adjustment
      if(!is.null(degree.bias) && degree.bias$method == 'SDS'){
        if(!multiplex && !heterogeneous){
          tmp = rwr_score / sds
          rwr_score = tmp / sum(tmp)
        }else{
          c_l = get.type(names(rwr_score), 2)
          for(j in seq_along(degree.bias$component)){
            id = which(c_l == degree.bias$component[j])
            rwr_score[id] = (rwr_score[id] / sds[id]) / (sum(rwr_score[id] / sds[id]) / sum(rwr_score[id]))
          }
        }
      }

      # Aggregate RWR scores
      if(flag){
        tmp.val = stats::aggregate(rwr_score, by = list(V(ig)$agg_name), FUN = agg.fun)
        rwr_score = tmp.val$x
        names(rwr_score) = tmp.val$Group.1
      }

      # Maximum Scoring Subgraph Algorithm
      rwr_score <- rwr_score - quantile(rwr_score, filtering_rate) # if filtering rate close to 0, quantile takes the minimum of rwr_score
      nets[[i]] <- heinz(ag, rwr_score)

      n.v <- vcount(nets[[i]])
      mZ <- sum(V(nets[[i]])$Z) / n.v
      mCC = sum(core_cc(nets[[i]])) / n.v

      # To prevent filtering out too many nodes at one iteration
      # This is the difference between the observed and theoretical filtering rates
      if((vcount(ag) - n.v) / vcount(ag) - filtering_rate >= 0.25){
        rm.id = c(rm.id, i)
        next
      }

      score_metrics[i,] <- c(mCC, mZ, n.v, grid[i])
    }
  }
  if(in.parallel){
    cl = parallel::makeForkCluster(n.cores, outfile = "")
    doParallel::registerDoParallel(cl)
    res = foreach::foreach(i = 1:length(grid), .verbose = FALSE, .packages = c("igraph", "Matrix"), .export = NULL, .noexport = NULL) %dopar% {
      # RWR
      PTmatrix2 <- RWR(nadjM = n.adj.M, setSeeds = seeds, restart = grid[i], heterogeneous = heterogeneous, multiplex = multiplex, net.weight = net.weight, layer.weight = layer.weight)
      rwr_score <- as.vector(PTmatrix2)
      names(rwr_score) <- V(ig)$name

      # Degree bias adjustment
      if(!is.null(degree.bias) && degree.bias$method == 'SDS'){
        if(!multiplex && !heterogeneous){
          tmp = rwr_score / sds
          rwr_score = tmp / sum(tmp)
        }else{
          c_l = get.type(names(rwr_score), 2)
          for(j in seq_along(degree.bias$component)){
            id = which(c_l == degree.bias$component[j])
            rwr_score[id] = (rwr_score[id] / sds[id]) / (sum(rwr_score[id] / sds[id]) / sum(rwr_score[id]))
          }
        }
      }

      # Aggregate RWR scores
      if(flag){
        tmp.val = stats::aggregate(rwr_score, by = list(V(ig)$agg_name), FUN = agg.fun)
        rwr_score = tmp.val$x
        names(rwr_score) = tmp.val$Group.1
      }

      # Maximum Scoring Subgraph Algorithm
      rwr_score <- rwr_score - quantile(rwr_score, filtering_rate) # if filtering rate close to 0, quantile takes the minimum of rwr_score
      nets <- heinz(ag, rwr_score)

      n.v <- vcount(nets)
      mZ <- sum(V(nets)$Z) / n.v
      mCC = sum(core_cc(nets)) / n.v

      rm.id = NULL

      # To prevent filtering out too many nodes at one iteration
      # This is the difference between the observed and theoretical filtering rates
      if((vcount(ag) - n.v) / vcount(ag) - filtering_rate >= 0.25){
        rm.id = i
        n.v = mZ = mCC = 0
      }

      list(c(mCC, mZ, n.v, grid[i]), nets, rm.id)
    }
    parallel::stopCluster(cl)
    score_metrics = do.call(rbind, lapply(res, function(x) x[[1]]))
    nets = lapply(res, function(x) x[[2]])
    rm.id = unlist(lapply(res, function(x) x[[3]]))
  }

  rm.id = unique(c(rm.id, which(score_metrics[,1] %in% c(0, NA))))
  if(length(rm.id) != 0){
    score_metrics = score_metrics[-rm.id,]
    nets = nets[-rm.id]
  }

  if(isTRUE(nrow(score_metrics) >= 1)){ # if score matrix isn't empty
    max.id <- which.max(score_metrics[, 2] * score_metrics[, 1])
    new.scores <- score_metrics[max.id,]
    new.g <- nets[[max.id]]
    best.r <- score_metrics[max.id, 4]
  }else{
    new.g = ag
    n.v = vcount(new.g)
    new.scores = c(sum(core_cc(new.g)) / n.v, sum(V(new.g)$Z) / n.v, n.v, NA)
    best.r = NA
  }

  return(list(new.g, new.scores, best.r))
}

#' @title Set the starting filtering rate
#'
#' @description
#' Sets the starting filtering rate (denoted as eta, the initial quantile of RWR scores used to shift scores before use in `heinz()`) in `run_AMEND()`.
#'
#' eta is linearly increasing with D(N,n), where D is some distance metric between the input network size N and the final subnetwork size n.
#' For method = "difference", D(N,n) = N - n
#' For method = "ratio", D(N,n) = N / n
#' For method = "log.ratio", D(N,n) = log(N / n)
#' Below are the boundaries within which eta changes linearly with D(N,n). Set to capture "normal" or "expected" input network sizes and user-defined final module sizes.
#' N: (500, 15000), n: (5, 200), eta: (0.3, 0.8)
#'
#' @param g input graph
#' @param n Final module size to approximate
#' @param method One of 'difference', 'ratio', or 'log.ratio'
#'
#' @returns Starting filtering rate
#'
get.eta0 = function(g, n, method = c("difference", "ratio", "log.ratio")){ # Function for setting the starting filtering rate
  eta.l = 0.3
  eta.u = 0.8
  N.l = 500
  N.u = 15000
  n.l = 5
  n.u = 200
  if(method == "difference"){
    D.l = N.l - n.u
    D.u = N.u - n.l
    D = vcount(g) - n
  }else if(method == "ratio"){
    D.l = N.l / n.u
    D.u = N.u / n.l
    D = vcount(g) / n
  }else if(method == "log.ratio"){
    D.l = log(N.l / n.u)
    D.u = log(N.u / n.l)
    D = log(vcount(g) / n)
  }else stop("method should be one of \'difference\', \'ratio\', or \'log.ratio\'")
  # Getting slope and intercept
  a = (eta.u - eta.l) / (D.u - D.l)
  b = (eta.u * D.l - eta.l * D.u) / (D.l - D.u)

  a * D + b
}

#' @title Exponential decay schedule for filtering rate
#'
#' @description Exponential decay schedule for filtering rate. Used to shift the raw RWR scores
#'
#' @param eta0 Starting filtering rate
#' @param d Decay parameter
#'
#' @returns vector of filtering rates for each iteration of AMEND
#'
#' @examples
#' rates = AMEND:::exp_filtering_rate(eta0 = 0.5, d = 0.2)
#' plot(1:100, rates)
#'
exp_filtering_rate <- function(eta0 = 0.5, d = 0.1){
  iterations <- 1:100
  return( eta0*exp(-d*(iterations - 1)) )
}

#' @title Sets the decay value
#'
#' @description Finds the max decay value for an exponential filtering rate schedule that will
#' allow the algorithm to arrive at a graph of size n.
#'
#' @param N Size of graph
#' @param eta0 Starting filtering rate
#' @param n Approximate size of final module
#' @param rate.diff Vector of differences between observed filtering rate and theoretical filtering rate
#'
#' @returns decay value, numeric
#'
#' @examples
#' d = AMEND:::find_decay(N = 1000, eta0 = 0.5, n = 50)
#' d
#'
find_decay <- function(N, eta0 = 0.5, n = 25, rate.diff = 0.07){
  rate.diff = mean(rate.diff)

  iter <- 1:100
  d.grid <- seq(0.01, 1, 0.01)
  decay <- c()
  k.max = 0.1
  k.min = -0.1

  for(j in 1:length(d.grid)){
    under <- 0
    perc <- exp_filtering_rate(eta0 = eta0, d = d.grid[j])

    d1 = N
    for(i in iter){
      #= EDIT (1/19/23): This will correct behavior when perc[i] is close to 1.
      # K = min(k.max, rate.diff / perc[i]) # previous assignment for K
      # rate difference up, K up, AFR up, decay up, theoretical filtering rate down, larger networks
      # rate.diff = observed minus theoretical filtering rate, so more nodes filtered out than expected, so need to push on the brakes i.e., increase decay.
      K = max(min(k.max, (rate.diff * (1 - perc[i])) / perc[i]), k.min) # Will lead to slightly lower adjusted filtering rate (afr) --> lower decay --> larger filtering rates --> smaller networks
      afr = perc[i] * (1 + K) # Adjusted filtering rate (AFR)
      d2 = round(d1 * (1 - afr), 0)
      if(d2 <= n){
        under <- 1
        break
      }
      if(d1 == d2){
        break
      }
      d1 = d2
    }
    if(under){
      decay <- c(decay, d.grid[j])
    }else break
  }

  if(length(decay) == 0){
    return(min(d.grid))
  }else return(max(decay))
}

#' @title Calculate the core-clustering coefficients of nodes in a graph
#'
#' @description `core_cc()` calculates the core-clustering coefficients of all nodes in a graph. This network concept was introduced by Bader and Hogue (2003)
#'
#' @details The core-clustering coefficient of a node is the density of the maximum k-core of the immediate neighborhood of that node. The k-core of a graph is the maximal subgraph in which every node has degree >= k.
#'
#' @param g Input graph. Must be connected
#' @param weight Logical. If true, the edge density is calculated taking into account edge weights
#'
#' @returns Vector of core-clustering coefficients for each node
#'
#' @examples
#' core_clust_coefs = AMEND:::core_cc(glut4_graph)
#' head(core_clust_coefs)
#'
core_cc = function(g, weight = TRUE){
  if(!"weight" %in% igraph::edge_attr_names(g) && weight){
    weight = FALSE
  }
  n = vcount(g)
  w = numeric(n)
  for(i in seq_len(n)){
    nborhood = igraph::make_ego_graph(g, order = 1, nodes = i, mode = "all")[[1]]
    cores = igraph::coreness(nborhood)
    max_core = igraph::induced_subgraph(nborhood, which(cores == max(cores)))
    w[i] = edge_density_weighted(max_core, weight = weight)
  }
  return(w)
}

#' @title Calculate the edge density of a weighted graph
#'
#' @description `edge_density_weighted()` calculates the edge density of a graph, taking into account edge weights.
#'
#' @details If weight is true, the edge density equals the sum of the edge weights divided by the maximum number of edges possible for the given graph. If weight is false, the edge density is the number of edges divided by the maximum number of edges possible for the given graph.
#'
#' @param graph Input graph
#' @param weight Logical. If true, the edge density is calculated taking into account edge weights
#'
#' @returns weighted edge density of graph
#'
#' @examples
#' AMEND:::edge_density_weighted(graph = glut4_graph, weight = TRUE)
#'
#' AMEND:::edge_density_weighted(graph = glut4_graph, weight = FALSE)
#'
#' # Compare to igraph::edge_density
#' igraph::edge_density(glut4_graph)
#'
#' # Remove weight attribute
#' g = igraph::delete_edge_attr(glut4_graph, "weight")
#' AMEND:::edge_density_weighted(graph = g, weight = TRUE)
#'
edge_density_weighted = function(graph, weight){
  if(weight){
    ew = igraph::E(graph)$weight
    if(any(ew > 1)){
      warning("Some edge weights are greater than one. Computing unweighted edge density.")
      return(igraph::edge_density(graph))
    }
    N = igraph::vcount(graph)
    return(sum(ew) / (N * (N - 1) / 2))
  }else return(igraph::edge_density(graph))
}

#' @title Get a larger subnetwork from previous iterations of AMEND
#'
#' @description
#'
#' This function allows the user to examine intermediate subnetworks that were generated during 'run_AMEND()', since the algorithm may overshoot the approximate final module size, or the user may be interested in the parent subnetworks of the final module.
#'
#' @param ig Input graph
#' @param amend_object Result from `run_AMEND()`
#' @param k Iteration associated with the subnetwork you want to retrieve
#'
#' @returns igraph object of the subnetwork associated with the given iteration.
#'
#' @examples
#' # Attach igraph library
#' library(igraph)
#'
#' # Inspect the igraph object included in AMEND package.
#' # One can see it has vertex attributes name, symbol, and ECI
#' glut4_graph
#' head(V(glut4_graph)$ECI)
#'
#' # Use run_AMEND() with an igraph object with a vertex attribute matching data.type arg
#' subnet1 = run_AMEND(graph = glut4_graph, n = 25, data = "ECI", FUN = "shift_scale",
#'                     FUN.params = list(DOI = -1, w = 0.5), normalize = "degree")
#'
#' # Inspect the statistics of the intermediate subnetworks
#' subnet1$stats
#'
#' # Retrieve the subnetwork of iteration 7
#' subnet7 = get_subnetwork(ig = glut4_graph, amend_object = subnet1, k = 7)
#' subnet7
#'
#' @export
get_subnetwork = function(ig, amend_object, k){
  nodes = amend_object$subnetworks[[k]]
  igraph::induced_subgraph(ig, which(V(ig)$name %in% nodes))
}

#' @title Get a pre-defined aggregation function
#'
#' @description Given a character scalar, this returns a function that will aggregate values of common nodes.
#'
#' @param x character scalar
#'
#' @return aggregation function
#'
get_aggregate_method = function(x){
  if(x == "mean"){
    agg.fun = mean
  }else if(x == "median"){
    agg.fun = stats::median
  }else if(x == "sum"){
    agg.fun = sum
  }else if(x == "gmean"){
    agg.fun = function(y) (prod(y))^(1/length(y))
  }else if(x == "hmean"){
    agg.fun = function(y) length(y) / (sum(1/y))
  }
  return(agg.fun)
}

#' @title Create an aggregated graph
#'
#' @description
#' An aggregated graph, in a multiplex context, is when the layers of a multiplex component are collapsed such that there is only one layer. A layer is designated as the primary, and edges in this layer take priority over edges of other layers. The collapsed multiplex will contain all of the original nodes.
#'
#' The purpose is to aggregate a multiplex component before the subnetwork identification step (using `heinz()`) so that the algorithm considers the edges of the primary layer. This process involves aggregating node attribute values and fusing nodes together such that the fused node is adjacent to all of the edges of its predecessors.
#'
#' @param graph igraph object. The multiplex graph to be aggregated.
#' @param control A named list. The element 'primary' contains the name of the primary layer for a multiplex component to be used during subnetwork identification. The element 'agg.method' contains a character scalar referring to an aggregation function.
#'
#' @return igraph object
#'
#' @seealso [melt_graph()], [expand_graph()], [run_AMEND()]
#'
create_aggregated_graph = function(graph, control){
  # Identify the multiplex components to collapse
  comps = extract_string(control$primary, "_", 1)
  # For each of these components, find primary layer
  layers = extract_string(control$primary, "_", 2)

  # !!! Be careful about lack of agreement between name and node_type!
  # node names will reflect the new labeling (just component, no layer info), while node_type will keep original '|component_layer' labeling.

  for(i in seq_along(comps)){
    # Strip layer info from node labels in primary layer
    V(graph)$name[V(graph)$node_type == control$primary[i]] = paste(get.type(V(graph)$name[V(graph)$node_type == control$primary[i]], 1), comps[i], sep = "|")
    v.target = get.type(V(graph)$name[V(graph)$node_type == control$primary[i]], 1)
    other.layers = stats::na.omit(unique(get.type(V(graph)$name, 2)[get.type(V(graph)$name, 3) == comps[i] & get.type(V(graph)$name, 4) != layers[i]]))
    for(j in seq_along(other.layers)){
      # Identify primary and non-primary nodes
      ## primary
      pn = which(get.type(V(graph)$name, 1) %in% v.target & get.type(V(graph)$name, 2) == other.layers[j])
      v.pn = V(graph)$name[pn]
      ## non-primary
      npn = which(!get.type(V(graph)$name, 1) %in% v.target & get.type(V(graph)$name, 2) == other.layers[j])
      v.npn = V(graph)$name[npn]
      ## Identify and delete primary-primary edges
      # Identify
      el = igraph::as_edgelist(graph)
      pp.id = which(apply(el, 1, function(x) all(x %in% v.pn)))
      # Delete
      graph = igraph::delete_edges(graph = graph, edges = pp.id)
      # Strip layer information from node labels
      V(graph)$name[get.type(V(graph)$name, 2) == other.layers[j]] = paste(get.type(V(graph)$name[get.type(V(graph)$name, 2) == other.layers[j]], 1), comps[i], sep = "|")
    }
    # Melt graph, i.e., fuse common nodes such that only one of them becomes adjacent to the edges of the others, and then remove duplicate nodes
    graph = melt_graph(g = graph, agg.method = control$agg.method)
  }
  return(graph)
}

#' @title Expand a graph
#'
#' @description
#' Given an aggregated graph, returns its expanded counterpart, induced from the original, full graph.
#' Takes nodes from aggregated graph, identifies nodes that were in a collapsed multiplex component, duplicates them as necessary, and appends appropriate '|component_layer' information.
#' Returns an induced subgraph of ig that contains only node names matching nodes of aggregated graph.
#'
#' @param ig Input graph
#' @param ag Aggregated graph
#' @param control Named list. Information on aggregated graph.
#'
#' @return igraph object
#'
expand_graph = function(ig, ag, control){
  if(!is.null(control)){
    ag.names = V(ag)$name
    # is.ag = get.type(ag.names, 3) %in% extract_string(control$primary, "_", 1)
    # ag.names[is.ag] = get.type(ag.names[is.ag], 1)

    ig.names = V(ig)$name
    is.ig = get.type(ig.names, 3) %in% extract_string(control$primary, "_", 1)
    # ig.names[is.ig] = get.type(ig.names[is.ig], 1)
    ig.names[is.ig] = paste(get.type(ig.names[is.ig], 1), get.type(ig.names[is.ig], 3), sep='|')

    graph = igraph::induced_subgraph(graph = ig, vids = which(ig.names %in% ag.names))
  }else graph = ag
  graph
}


#' @title Set the Biased Random Walk attribute
#'
#' @description
#' Calculates values for a biased random walk. This value is a decreasing function of minimum distance in the full graph to a node of interest (NOI): f(d,k) = exp(-k*d), where k >= 0, d=distance.
#'
#' @param graph Input graph
#' @param brw.attr Character vector. Denotes the nodes of interest (NOI) used to artificially seed the active module.
#' @param brw.attr.nm Character scalar. The name of the vertex attribute to be used for biased random walk.
#' @param brw.flag Logical. True returns graph with new biased random walk attribute. False returns input graph.
#'
#' @return igraph object
#'
set_brw_attr = function(graph, brw.attr, brw.attr.nm, brw.flag){
  if(brw.flag){
    # c.id = which(unlist(lapply(brw.attr, is.character)))
    # n.id = which(unlist(lapply(brw.attr, is.numeric)))
    # noi = brw.attr[[c.id]]
    # k = brw.attr[[n.id]]
    noi = brw.attr
    k = 1
    # Get IDs of NOIs in graph
    name.only = get.type(V(graph)$name, 1)
    name.comp = paste(name.only, get.type(V(graph)$name, 3), sep = "|")
    name.comp.layer = V(graph)$name
    noi.id = which(name.only %in% noi | name.comp %in% noi | name.comp.layer %in% noi)
    if(length(noi.id) > 0){
      # Get distance to closest NOI for each node in graph
      d.mat = igraph::distances(graph = graph, v = noi.id)
      node.dists = apply(d.mat, 2, mean)
      # Calculate values for biased random walk
      alpha = exp(-node.dists * k)
      while(any(alpha == 0 & !node.dists %in% c(0, Inf))){
        k = k * 0.9
        alpha = exp(-node.dists * k)
      }
      # Assign values to 'brw.attr.nm' vertex attribute in graph
      igraph::vertex_attr(graph, brw.attr.nm) = alpha[match(V(graph)$name, names(alpha))]
    }
  }
  graph
}
