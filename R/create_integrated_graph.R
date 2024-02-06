#' @importFrom igraph vcount V V<- E E<- vertex_attr vertex_attr<-

#' @title Create an integrated graph from a list of graph-like objects
#'
#' @description
#' Given a single graph-like object (igraph, adjacency matrix, or edge list) or a list of these, create an integrated graph with seed values computed from _data_, _FUN_, and _FUN.params_.
#'
#' @details
#' Heterogeneous graphs are graphs with more than one node type. We call these different subnetworks of node types _components_. A multiplex graph has more than one edge or data type within a component, with the different subnetworks called _layers_.
#'
#' In order to be able to properly identify nodes in the network, there are a few naming conventions for both node labels and function arguments that must be followed.
#'
#' If the input network is heterogeneous and a single object, node names must be appended with '|component', where 'component' is the name of the component and '|' is the delimiter (so node names can't contain '|'). If a list input, the list names correspond to the component and will be appended to node names in the integrated graph. The names of list elements corresponding to bipartite networks must include the relevant components separated by ';', e.g., 'protein;metabolite;drug' (so names of components can't include ';').
#'
#' If the input network is multiplex and a single object, node names must be appended with '|component_layer', where 'layer' is the name of the layer and '`_`' is the delimiter (so component and layer names can't contain '_'). If a list input, the list names correspond either to the component _and_ layer, or if it only gives the component name and that component is multiplex, the node names must be appended with 'component_layer' to differentiate between layers. Bipartite connections between nodes of two multiplex components will be applied to all layers.
#'
#' An effort was made to flexibly accept different network input objects. Just keep in mind that for heterogeneous/multiplex networks, the algorithm must have a way to identify a specific node within a specific component and layer using the '|' and '_' syntax described above.
#'
#' @param graph,adj_matrix,edge_list A single graph like object in the form of an igraph, adjacency matrix, or edge list. Or, a named list containing multiple graph-like objects of the same type, to be merged. The merged graph must be connected. If a third column is provided in an edge list, these are taken as edge weights. Only one of graph, adj_matrix, or edge_list should be given, with priority given to graph, then adj_matrix, then edge_list. See 'Details' for correct naming conventions.
#' @param data A named list of named numeric vectors (list elements correspond to graph components), a named numeric vector, or a character scalar (denoting a vertex attribute of the input igraph object) containing the experimental data from which seed values for RWR will be derived according to _FUN_ and _FUN.params_ args. See 'Details' for correct naming conventions.
#' @param node_type A named list of character vectors (list elements correspond to graph components), a named character vector, a character scalar (denoting a vertex attribute of the input igraph object), or NULL. The list elements represent sets of nodes that belong to the corresponding components. The character vector contains _component_layer_ tags for the node name given in names of vector. If NULL and multiplex and/or heterogeneous, node labels must follow 'name|type_layer' naming scheme (e.g., MYC|gene_1). See 'Details' for correct naming conventions.
#' @param brw.attr A named list of named numeric vectors (list elements correspond to graph components), a named numeric vector, a character scalar (denoting a vertex attribute of the input igraph object), or NULL. Biased random walk vertex attribute values. Should be non-negative, with larger values increasing the transition probabilities to a node in RWR. If NULL, all nodes are given a value of 1. See 'Details' for biased random walk info.
#' @param FUN A function, named list of functions, named list of character scalars, a single character scalar, or NULL. Function for transforming values in _data_ to derive seed values for RWR. Names correspond to graph components. Character strings should correspond to a default function: one of 'binary', 'shift_scale', 'p_value', or 'exp'. NULL means no transformation is done to values in _data_. See 'Details' for descriptions of default functions.
#' @param FUN.params A named list of lists of named function arguments, a named list of named function arguments, or NULL. Function arguments to be passed to _FUN_. Names should match names in _FUN_.
#' @param heterogeneous Logical. If TRUE, graph is considered heterogeneous (more than one distinct node type, e.g., proteins and metabolites), and node_type must be included as an argument or graph vertex attribute.
#' @param multiplex Logical. If true, graph is assumed to contain multiplex components.
#' @param lcc Logical. If true, return Largest Connected Component
#'
#' @returns an igraph of the integrated graph with vertex attributes 'seeds' and 'Z' computed from data, FUN, and FUN.params
#'
#' @seealso [run_AMEND()], [list_input_checks()], [transition_matrix()], [RWR()], [heinz()], [get_subnetwork()]
#'
#' @examples
#' # example code
#' \dontrun{
#' #=== Monoplex-heterogeneous network ===#
#' # using a list of edge lists, with data as list input
#' # Create object for 'data' argument
#' data.list.input = vector("list", 2)
#' names(data.list.input) = c("prot", "meta")
#' for(i in seq_along(data.list.input)){
#'   uniq.names = extract_string(unique(c(edgelists[[i]][,1],
#'                                      edgelists[[i]][,2])), "\\|", 1)
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
#' g1 = create_integrated_graph(edge_list = edgelists, data = data.list.input,
#'                              node_type = NULL, brw.attr = NULL,
#'                              FUN = FUN, FUN.params = FUN.params,
#'                              heterogeneous = TRUE, multiplex = FALSE)
#'
#' #=== Multiplex-heterogeneous network with list of graphs, data as vertex attribute ==#
#' # Create data vertex attributes, called 'scores'
#' v.attr.name = "scores"
#' node.types = c("prot_1", "prot_2", "prot_3", "meta")
#' for(i in seq_along(node.types)){
#'   # p-value
#'   igraph::vertex_attr(list.of.graphs[[node.types[i]]], v.attr.name) =
#'   runif(vcount(list.of.graphs[[node.types[i]]]))
#' }
#'
#' # Setting multiplex/heterogeneous parameters
#' layer.weight = list(prot = rep(1/3, 3))
#' switch.layer.prob = list(prot = c(prot_1 = 0.5, prot_2 = 0.2, prot_3 = 0.8))
#' net.weight = c(prot = 0.5, meta = 0.5)
#' jump.prob = c(prot = 0.5, meta = 0.3)
#'
#' g2 = create_integrated_graph(graph = list.of.graphs, data = "scores", node_type = NULL,
#'                              brw.attr = NULL, FUN = "p_value", FUN.params = NULL,
#'                              heterogeneous = TRUE, multiplex = TRUE)
#'
#' #=== Multiplex-homogeneous network ===#
#' # list of adjacency matrices, data as list input
#' # Create a list of binary vectors
#' # (e.g., significantly up- or down-regulated genes)
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
#' g3 = create_integrated_graph(adj_matrix = list.of.adjmats, data = data.list.input,
#'                              node_type = NULL, brw.attr = NULL, FUN = "binary",
#'                              FUN.params = NULL, heterogeneous = FALSE, multiplex = TRUE)
#' }
#' @export
create_integrated_graph = function(graph = NULL, adj_matrix = NULL, edge_list = NULL, data, node_type = NULL, brw.attr = NULL,
                                   FUN = NULL, FUN.params = NULL, heterogeneous = FALSE, multiplex = FALSE, lcc = FALSE){
  # Handling data and brw.attr arguments
  if(is.character(data) && length(data) == 1){
    data.name = data
  }else data.name = "pre.seed.values"
  if(is.character(brw.attr) && length(brw.attr) == 1){
    brw.attr.nm = brw.attr
  }else brw.attr.nm = "brw.values"

  # The goal is to get all node name, node type/layer, and edge information into an igraph object for easy argument checking and easy switch to adjacency mat,
  # with functions able to manipulate these adj mats for transition matrix and RWR.
  # Also quickly need to establish Unique node types, which are mono or multiplex, layers of each node type, easy functions for quick access to this info based on 'name|type_layer' naming scheme.
  # To merge all edge and node network information (not v.attrs), I convert to an edgelist all components, then I rbind if necessary and convert to an igraph object, removing loops and multiple edges.
  # From here I check connectivity, do name checks, then v.attr & function-argument checks
  # If multiplex, convert adjacency matrix to connect layers.

  # Handling non-graph inputs
  if(is.null(graph)){
    if(!is.null(adj_matrix)){
      if(is.list(adj_matrix) && !is.null(names(adj_matrix))){
        # Get a named list of node names in each list element
        adj_matrix = lapply(adj_matrix, adj_mat_checks)
        node_names = lapply(adj_matrix, rownames)
        node_names = list_input_checks(nn = node_names, multiplex = multiplex, heterogeneous = heterogeneous)
        for(i in seq_along(node_names)){
          rownames(adj_matrix[[i]]) = colnames(adj_matrix[[i]]) = node_names[[i]]
        }
        edge_list = lapply(adj_matrix, function(x){
          x = igraph::graph_from_adjacency_matrix(adjmatrix = x, mode = "undirected", weighted = TRUE)
          el = igraph::as_edgelist(graph = x, names = TRUE)
          if("weight" %in% igraph::edge_attr_names(x)){
            el =  cbind(el, igraph::E(x)$weight)
          }else el = cbind(el, 1)
          el
        })
        for(i in seq_along(edge_list)){
          if(i == 1){
            el = edge_list[[i]]
          }else{
            el = rbind(el, edge_list[[i]])
          }
        }
        graph = igraph::graph_from_edgelist(el = el[,1:2], directed = FALSE)
        igraph::E(graph)$weight = as.numeric(el[,3])
        graph = igraph::simplify(graph = graph, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "median")
        # rm(adj_matrix, edge_list, el, node_names)
      }else{
        adj_matrix = lapply(list(adj_matrix), adj_mat_checks)[[1]]
        graph = igraph::graph_from_adjacency_matrix(adjmatrix = adj_matrix, mode = "undirected", weighted = TRUE)
        # rm(adj_matrix)
      }
    }else if(!is.null(edge_list)){
      if(is.list(edge_list) && !is.null(names(edge_list))){
        # Get a named list of node names in each list element
        edge_list = lapply(edge_list, el_checks)
        graph = lapply(edge_list, function(x){
          g = igraph::graph_from_edgelist(el = x[,1:2], directed = FALSE)
          igraph::E(g)$weight = as.numeric(x[,3])
          g = igraph::simplify(graph = g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "median")
          g
        })
        node_names = lapply(graph, function(x) V(x)$name)
        node_names = list_input_checks(nn = node_names, multiplex = multiplex, heterogeneous = heterogeneous)
        for(i in seq_along(node_names)){
          V(graph[[i]])$name = node_names[[i]]
        }
        edge_list = lapply(graph, function(x){
          el = igraph::as_edgelist(graph = x, names = TRUE)
          if("weight" %in% igraph::edge_attr_names(x)){
            el =  cbind(el, igraph::E(x)$weight)
          }else el = cbind(el, 1)
          el
        })
        for(i in seq_along(edge_list)){
          if(i == 1){
            el = edge_list[[i]]
          }else{
            el = rbind(el, edge_list[[i]])
          }
        }
        graph = igraph::graph_from_edgelist(el = el[,1:2], directed = FALSE)
        igraph::E(graph)$weight = as.numeric(el[,3])
        graph = igraph::simplify(graph = graph, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "median")
        # rm(adj_matrix, edge_list, el, node_names)
      }else{
        edge_list = lapply(list(edge_list), el_checks)[[1]]
        graph = igraph::graph_from_edgelist(el = edge_list[,1:2], directed = FALSE)
        igraph::E(graph)$weight = as.numeric(edge_list[,3])
        graph = igraph::simplify(graph = graph, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "median")
        # rm(edge_list)
      }
    }
  }else if(is.list(graph) && !is.null(names(graph))){
    # Handling list of igraph objects
    # node_type v.attr takes priority over node_type info in name or in node_type arg if provided.
    # If there is a v.attr 'data.name', extract as a named vector. Assumes that 'data' was a character string
    # If there is a v.attr 'brw.attr.nm', extract this as a named vector. Assumes that 'brw.attr' was a character string
    # node_type is a named list of character vectors. Each entry name corresponds to a distinct node type & layer (if one exists)
    # data can be a named list of named numeric vectors. Entry names correspond to node types, and optionally may refer to layers (if 1 layer given, all must be)(If none given and multiplex, map same values to all layers).
    # brw.attr is analogous to data
    graph = lapply(graph, graph_checks)
    node_names = lapply(graph, function(x) V(x)$name)
    node_names = list_input_checks(nn = node_names, multiplex = multiplex, heterogeneous = heterogeneous)
    for(i in seq_along(graph)) V(graph[[i]])$name = node_names[[i]]

    ## Collect any vertex attributes from graphs
    v.attrs = unique(unlist(lapply(graph, function(x) igraph::vertex_attr_names(x))))
    if(is.character(data) && length(data) == 1){
      v.attrs = unique(c(v.attrs, data.name))
    }
    if(is.character(brw.attr) && length(brw.attr) == 1){
      v.attrs = unique(c(v.attrs, brw.attr))
    }
    if(length(v.attrs) != 0){
      tmp.scores = vector("list", length(v.attrs)); names(tmp.scores) = v.attrs
      for(i in seq_along(v.attrs)){
        tmp.scores[[i]] = unlist(lapply(unname(graph), function(x){
          if(v.attrs[i] %in% igraph::vertex_attr_names(x)){
            a = igraph::vertex_attr(x, v.attrs[i])
            names(a) = V(x)$name
            a
          }else c()
        }))
        tmp.scores[[i]] = tmp.scores[[i]][match(unique(names(tmp.scores[[i]])), names(tmp.scores[[i]]))]
      }
    }

    edge_list = lapply(graph, function(x){
      el = igraph::as_edgelist(graph = x, names = TRUE)
      if("weight" %in% igraph::edge_attr_names(x)){
        el =  cbind(el, igraph::E(x)$weight)
      }else el = cbind(el, 1)
      el
    })
    for(i in seq_along(edge_list)){
      if(i == 1){
        el = edge_list[[i]]
      }else{
        el = rbind(el, edge_list[[i]])
      }
    }
    graph = igraph::graph_from_edgelist(el = el[,1:2], directed = FALSE)
    igraph::E(graph)$weight = as.numeric(el[,3])
    graph = igraph::simplify(graph = graph, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "median")
    if(length(v.attrs) != 0){
      # Redistribute the vertex attributes to the integrated graph
      for(i in seq_along(v.attrs)){
        ind = match(names(tmp.scores[[i]]), V(graph)$name)
        igraph::vertex_attr(graph, v.attrs[i], ind[!is.na(ind)]) = unname(tmp.scores[[i]][!is.na(ind)])
        # Switch NAs to zero
        if(is.numeric(igraph::vertex_attr(graph, v.attrs[i]))) igraph::vertex_attr(graph, v.attrs[i], which(is.na(igraph::vertex_attr(graph, v.attrs[i])))) = 0
      }
    }
  }else if(!igraph::is_igraph(graph)) stop("Unrecognized input type for 'graph', 'adj_matrix', or 'edge_list'.")

  #==================================#
  # BEGIN single igraph object checks
  #==================================#

  if(!"name" %in% igraph::vertex_attr_names(graph)) stop("No 'name' vertex attribute.")

  # Check node_type argument and assign unique IDs as 'name|type_layer' or 'name|type' for monoplex components.
  no.type.id = which(!grepl("\\|", V(graph)$name))
  if(length(no.type.id) == 0){ # If they all have '|', which must be the case in a multiplex-homogeneous context, but this being True doesn't require it to be multi-homo. list_input_checks() should verify that proper naming conventions were followed, as the rest of the function depends on this!
    if(any(n.char(V(graph)$name, "\\|") != 1)) stop("Incorrect use of \"|\" in node labels. See function documentation for details on naming conventions.")
    if(!"node_type" %in% igraph::vertex_attr_names(graph)) V(graph)$node_type = get.type(V(graph)$name,2)
  }else{
    # This 'else' section doesn't include mulit-homo graphs since this will be captured by if(length(no.type.id) == 0)
    if(!"node_type" %in% igraph::vertex_attr_names(graph)){
      # Assuming that all multiplex nodes (except those representing bipartite edges) have 'name|type_layer' format.
      # So in multi-hetero context, the elements in 'no.type.id' will represent nodes of bipartite edges whose node type is multiplex, or monoplex nodes in hetero/homo graph
      if(is.null(node_type) && heterogeneous) stop("There is no vertex attribute \'node_type\' or node_type argument. At least one must be given for heterogeneous graphs.")
      if(heterogeneous && !multiplex){ # If mono-hetero
        if(is.character(node_type) && !is.null(names(node_type))){ # If a named character vector
          if(any(!V(graph)$name[no.type.id] %in% names(node_type))) stop("node_type does not contain all nodes in graph.")
          # Append component to node names ('name|component_layer')
          V(graph)$name[no.type.id] = paste(V(graph)$name[no.type.id], node_type[match(V(graph)$name[no.type.id], names(node_type))], sep = "|")
        }else if(is.list(node_type) && !is.null(names(node_type))){ # If named list
          nm = V(graph)$name[no.type.id]
          tmp = character(length(no.type.id))
          for(i in seq_along(nm)){
            t = names(node_type)[unlist(lapply(node_type, function(x) nm[i] %in% x))]
            if(length(t) == 0) t = "" # stop("node_type does not contain all nodes in graph.")
            tmp[i] = t
          }
          V(graph)$name[no.type.id] = paste(nm, tmp, sep = "|")
        }else stop("Unrecognized input for 'node_type'.")
      }else if(heterogeneous && multiplex){
        if(is.list(node_type) && !is.null(names(node_type))){ # If named list
          nm = V(graph)$name[no.type.id]
          tmp = character(length(no.type.id))
          for(i in seq_along(nm)){
            t = names(node_type)[unlist(lapply(node_type, function(x) nm[i] %in% x))]
            if(length(t) == 0) t = "_" # stop("node_type does not contain all nodes in graph.")
            if(length(t) > 1) if(!(all(n.char(t, "_") == 1) && length(unique(extract_string(t, "_", pos=1))) == 1 && length(unique(extract_string(t, "_", pos=2))) == length(t))) stop("Some nodes map to multiple node types.")
            t = extract_string(t, "_", pos=1)[1]
            tmp[i] = t
          }
          V(graph)$name[no.type.id] = paste(nm, tmp, sep = "|")
        }else stop("Unrecognized input for 'node_type'.")
      }else if(!heterogeneous && !multiplex){
        # Assign arbitrary 'name|type' label for monoplex-homogeneous graphs, for uniform handling of objects throughout the function
        V(graph)$name = paste(V(graph)$name, "typeA", sep = "|")
      }
      # Assign node_type vertex attribute
      V(graph)$node_type = get.type(V(graph)$name,2)
    }else V(graph)$name[no.type.id] = paste(V(graph)$name[no.type.id], V(graph)$node_type[no.type.id], sep = "|") # If there is a 'node_type' v.attr
  }

  # Remove any names with "" after "|", if heterogeneous or multiplex. These are any nodes that didn't map to a node type
  if(heterogeneous || multiplex) graph = igraph::delete_vertices(graph = graph, v = which(is.na(V(graph)$node_type)))

  # Melt the graph (i.e., convert to edge list and back to graph) to fuse nodes from bipartite graphs to their corresponding monoplex components. Conserve any v.attrs
  graph = melt_graph(g = graph)

  # Check uniqueness of all 'name|type_layer' node IDs
  if(length(unique(V(graph)$name)) < igraph::vcount(graph)) stop("Not all node names are unique.")

  uniq.types = unique(V(graph)$node_type)
  is.multi = any(grepl("_", uniq.types))
  is.hetero = length(unique(extract_string(uniq.types, "_", pos=1))) > 1

  if(multiplex && !is.multi) stop("multiplex=TRUE but no multiplex layers found. See function documentation for the correct node type naming conventions.")
  if(heterogeneous && !is.hetero) stop("heterogeneous=TRUE but only one unique node type found. See function documentation for the correct node type naming conventions.")
  if(!multiplex && is.multi) stop("multiplex=FALSE but multiplex layers were found. See function documentation for the correct node type naming conventions.")
  if(!heterogeneous && is.hetero) stop("heterogeneous=FALSE but more than one node type found. See function documentation for the correct node type naming conventions.")


  ### Multiplex or Heterogeneous
  if(multiplex || heterogeneous){
    if(multiplex){
      multiplex.comps = unique(extract_string(uniq.types[grepl("_", uniq.types)], "_", pos=1)) # NULL if no '_
      bp.nm = uniq.types[uniq.types %in% multiplex.comps & !grepl("_", uniq.types)]
    }else{
      bp.nm = NULL
    }
    non.bp.nm = setdiff(uniq.types, bp.nm) # node types that possibly map to elements in function args

    # 'data'
    if(is.character(data) && length(data) == 1){
      if(!data.name %in% igraph::vertex_attr_names(graph)) stop(paste0("No vertex attribute '", data.name,"' found."))
    }else if(is.list(data) && !is.null(names(data))){ # Map list values to graph
      for(i in seq_along(non.bp.nm)){
        if(!non.bp.nm[i] %in% names(data) && !extract_string(non.bp.nm[i], "_", pos=1) %in% names(data)) stop("Node types don't match with names of 'data'.")
        if(grepl("_", non.bp.nm[i]) && non.bp.nm[i] %in% names(data) && extract_string(non.bp.nm[i], "_", pos=1) %in% names(data)) stop("Incorrect labeling for 'data'. See function documentation for correct naming conventions.")
        t = V(graph)$node_type == non.bp.nm[i]
        if(non.bp.nm[i] %in% names(data)){
          ind = match(extract_string(V(graph)$name, "\\|", pos=1), names(data[[match(non.bp.nm[i], names(data))]]))
          igraph::vertex_attr(graph, data.name, which(!is.na(ind))) = unname(data[[match(non.bp.nm[i], names(data))]][ind[!is.na(ind)]])
          if(any(is.na(ind) & V(graph)$node_type %in% c(non.bp.nm[i], extract_string(non.bp.nm[i], "_", pos=1)))) warning(paste0(sum(is.na(ind) & V(graph)$node_type %in% c(non.bp.nm[i], extract_string(non.bp.nm[i], "_", pos=1))), " of node type ", non.bp.nm[i], " didn't map to 'data'. These will be given seed values of 0."))
        }else if(!non.bp.nm[i] %in% names(data) && extract_string(non.bp.nm[i], "_", pos=1) %in% names(data)){
          ind = match(extract_string(V(graph)$name, "\\|", pos=1), names(data[[match(extract_string(non.bp.nm[i], "_", pos=1), names(data))]]))
          ind[V(graph)$node_type != non.bp.nm[i]] = NA
          igraph::vertex_attr(graph, data.name, which(!is.na(ind))) = unname(data[[which(names(data) == extract_string(non.bp.nm[i], "_", pos=1))]][ind[!is.na(ind)]])
          if(any(is.na(ind[V(graph)$node_type == non.bp.nm[i]]) & V(graph)$node_type[V(graph)$node_type == non.bp.nm[i]] %in% c(non.bp.nm[i], extract_string(non.bp.nm[i], "_", pos=1)))) warning(paste0(sum(is.na(ind[V(graph)$node_type == non.bp.nm[i]]) & V(graph)$node_type[V(graph)$node_type == non.bp.nm[i]] %in% c(non.bp.nm[i], extract_string(non.bp.nm[i], "_", pos=1))), " of node type ", non.bp.nm[i], " didn't map to 'data'. These will be given seed values of 0."))
        }else stop("Incorrect labeling for 'data'. See function documentation for correct naming conventions.")
      }
    }else if(is.numeric(data) && !is.null(names(data))){ # Map same values to same node across different layers. NAs for nodes that don't match (turn to 0s for seed & Z)
      ind = match(extract_string(V(graph)$name, "\\|", pos=1), names(data))
      if(any(is.na(ind))) warning(paste0(sum(is.na(ind)), " nodes didn't map to 'data'. These will be given seed values of 0."))
      igraph::vertex_attr(graph, data.name, which(!is.na(ind))) = unname(data[ind[!is.na(ind)]])
    }else stop("Incorrect input for 'data'.")

    # 'brw.attr'
    if(is.character(brw.attr) && length(brw.attr) == 1){
      if(!brw.attr.nm %in% igraph::vertex_attr_names(graph)) stop(paste0("No vertex attribute '", brw.attr.nm,"' found."))
    }else if(is.list(brw.attr) && !is.null(names(brw.attr)) && !(length(brw.attr) == 2 && any(unlist(lapply(brw.attr, is.character))) && any(unlist(lapply(brw.attr, is.numeric))))){ # Map list values to graph
      for(i in seq_along(non.bp.nm)){
        if(!non.bp.nm[i] %in% names(brw.attr) && !extract_string(non.bp.nm[i], "_", pos=1) %in% names(brw.attr)) stop("Node types don't match with names of 'brw.attr'.")
        if(grepl("_", non.bp.nm[i]) && non.bp.nm[i] %in% names(brw.attr) && extract_string(non.bp.nm[i], "_", pos=1) %in% names(brw.attr)) stop("Incorrect labeling for 'brw.attr'. See function documentation for correct naming conventions.")
        if(non.bp.nm[i] %in% names(brw.attr)){
          ind = match(extract_string(V(graph)$name, "\\|", pos=1), names(brw.attr[[match(non.bp.nm[i], names(brw.attr))]]))
          igraph::vertex_attr(graph, brw.attr.nm, which(!is.na(ind))) = unname(brw.attr[[match(non.bp.nm[i], names(brw.attr))]][ind[!is.na(ind)]])
          if(any(is.na(ind) & V(graph)$node_type %in% c(non.bp.nm[i], extract_string(non.bp.nm[i], "_", pos=1)))) warning(paste0(sum(is.na(ind) & V(graph)$node_type %in% c(non.bp.nm[i], extract_string(non.bp.nm[i], "_", pos=1))), " of node type ", non.bp.nm[i], " didn't map to 'brw.attr'. These will be given values of 1."))
        }else if(!non.bp.nm[i] %in% names(brw.attr) && extract_string(non.bp.nm[i], "_", pos=1) %in% names(brw.attr)){
          ind = match(extract_string(V(graph)$name, "\\|", pos=1), names(brw.attr[[match(extract_string(non.bp.nm[i], "_", pos=1), names(brw.attr))]]))
          ind[V(graph)$node_type != non.bp.nm[i]] = NA
          igraph::vertex_attr(graph, brw.attr.nm, which(!is.na(ind))) = unname(brw.attr[[which(names(brw.attr) == extract_string(non.bp.nm[i], "_", pos=1))]][ind[!is.na(ind)]])
          if(any(is.na(ind[V(graph)$node_type == non.bp.nm[i]]) & V(graph)$node_type[V(graph)$node_type == non.bp.nm[i]] %in% c(non.bp.nm[i], extract_string(non.bp.nm[i], "_", pos=1)))) warning(paste0(sum(is.na(ind[V(graph)$node_type == non.bp.nm[i]]) & V(graph)$node_type[V(graph)$node_type == non.bp.nm[i]] %in% c(non.bp.nm[i], extract_string(non.bp.nm[i], "_", pos=1))), " of node type ", non.bp.nm[i], " didn't map to 'brw.attr'. These will be given values of 1."))
        }else stop("Incorrect labeling for 'brw.attr'. See function documentation for correct naming conventions.")
      }
      igraph::vertex_attr(graph, brw.attr.nm, which(is.na(igraph::vertex_attr(graph, brw.attr.nm)))) = 1
    }else if(is.numeric(brw.attr) && !is.null(names(brw.attr))){ # Map same values to same node across different layers. NAs for nodes that don't match (turn to 0s for seed & Z)
      ind = match(extract_string(V(graph)$name, "\\|", pos=1), names(brw.attr))
      if(any(is.na(ind))) warning(paste0(sum(is.na(ind)), " nodes didn't map to 'brw.attr'. These will be given values of 1 (unbiased RW)."))
      igraph::vertex_attr(graph, brw.attr.nm, which(!is.na(ind))) = brw.attr[ind[!is.na(ind)]]
      igraph::vertex_attr(graph, brw.attr.nm, which(is.na(igraph::vertex_attr(graph, brw.attr.nm)))) = 1
    }else if(is.null(brw.attr)){ # Default to 1 (equivalent to an unbiased RW)
      igraph::vertex_attr(graph, brw.attr.nm) = 1
    }else if(!(is.character(brw.attr) && length(brw.attr) > 1)) stop("Incorrect input for 'brw.attr'.")

    # Check FUN argument
    if(is.null(FUN) || (is.character(FUN) && length(FUN) == 1) || (is.function(FUN) && length(FUN) == 1)){
      # Assign actual function according to names...
      if(!is.function(FUN)) FUN = get_default_function(FUN)
      fun.map = FALSE
    }else if(is.list(FUN) && !is.null(names(FUN))){
      # Checking names
      fun.nm = unique(names(FUN))
      for(i in seq_along(fun.nm)) if(sum(names(FUN) == fun.nm[i]) > 1) stop("Non-unique names for 'FUN'.")
      FUN = lapply(FUN, function(x){
        if(is.character(x) && length(x) == 1){
          x = get_default_function(x)
        }else if(!is.function(x)){
          stop("Unrecognized input for 'FUN'.")
        }
        x
      })
      fun.map = TRUE
    }else stop("Incorrect input type for 'FUN'.")

    # Check FUN.params argument
    if(is.null(FUN.params)){
      fp.map = FALSE
    }else if(is.list(FUN.params) && !is.null(names(FUN.params))){
      if(all(unlist(lapply(FUN.params, function(x) !is.list(x))))){
        fp.map = FALSE
      }else if(all(unlist(lapply(FUN.params, function(x) is.list(x))))){
        # Checking names
        if(any(unlist(lapply(FUN.params, function(x) is.null(names(x)))))) stop("Elements in FUN.params should be lists of named elements.")
        fun.nm = unique(names(FUN.params))
        for(i in seq_along(fun.nm)) if(sum(names(FUN.params) == fun.nm[i]) > 1) stop("Names of 'graph' map to more than one name in 'FUN.params'.")
        fp.map = TRUE
      }else stop("Incorrect input type for 'FUN.params'.")
    }else stop("Incorrect input type for 'FUN.params'.")

    # For each non-BP component of graph, apply the corresponding function from FUN with corresponding args from FUN.params. Check function results
    seeds = Z = numeric(vcount(graph))
    for(i in seq_along(non.bp.nm)){
      ids = which(V(graph)$node_type == non.bp.nm[i])
      dat = igraph::vertex_attr(graph, data.name, ids)
      # Mapping graph name to FUN (if it is a named list of lists)
      if(fun.map){
        new.nm = ifelse(!non.bp.nm[i] %in% names(FUN), extract_string(non.bp.nm[i], "_", pos=1), non.bp.nm[i])
        id.f = which(names(FUN) %in% new.nm)
        if(length(id.f) == 0){ # If graph doesn't map to any name in FUN
          warning("graph doesn't map to any name of 'FUN'. Not applying any function.")
          use.f = function(x,...) x
        }else use.f = FUN[[id.f]]
      }else use.f = FUN
      # Mapping graph name to FUN.params (if it is a named list of lists or of function args)
      if(fp.map){
        new.nm = ifelse(!non.bp.nm[i] %in% names(FUN.params), extract_string(non.bp.nm[i], "_", pos=1), non.bp.nm[i])
        id.fp = which(names(FUN.params) %in% new.nm)
        if(length(id.fp) == 0){ # If graph doesn't map to any name in FUN.params
          use.fp = NULL
        }else use.fp = FUN.params[[id.fp]]
      }else use.fp = FUN.params
      # Only pass function args if non-NULL
      if(is.null(use.fp)){
        seeds[ids] = do.call(use.f, list(dat))
      }else seeds[ids] = do.call(use.f, c(list(dat), use.fp))
      seeds[ids][is.na(seeds[ids])] = 0
      Z[ids] = (seeds[ids] - mean(seeds[ids])) / sd(seeds[ids])
      if(!all.good(seeds[ids])) stop("Applying 'FUN' to 'data' resulted in either negative or all zero values.")
    }
    igraph::vertex_attr(graph, "seeds") = seeds
    igraph::vertex_attr(graph, "Z") = Z

  }else{ ### Monoplex-homogeneous
    v.nm = V(graph)$name
    # Check data argument
    if(is.character(data) && length(data) == 1){
      if(!data %in% igraph::vertex_attr_names(graph)) stop(paste0("No vertex attribute \'", data, "\' in graph."))
    }else if(is.numeric(data) && !is.null(names(data))){
      if(any(!extract_string(v.nm, "\\|", pos=1) %in% names(data))) warning(paste0(sum(!extract_string(v.nm, "\\|", pos=1) %in% names(data)), " nodes didn't map to 'data'. Unmapped nodes will be given seed values of 0."))
      ind = match(names(data), extract_string(v.nm, "\\|", pos=1))
      vertex_attr(graph, data.name, ind[!is.na(ind)]) = data[!is.na(ind)]
    }else stop("Incorrect input for 'data'.")
    # Check brw.attr argument
    if(is.null(brw.attr)){
      igraph::vertex_attr(graph, brw.attr.nm) = 1
    }else if(is.character(brw.attr) && length(brw.attr) == 1){
      if(!brw.attr %in% igraph::vertex_attr_names(graph)){
        warning(paste0("No vertex attribute \'", brw.attr, "\' in graph. Continuing with unbiased RW."))
        igraph::vertex_attr(graph, brw.attr.nm) = 1
      }
    }else if(is.numeric(brw.attr) && !is.null(names(brw.attr))){
      if(any(!extract_string(v.nm, "\\|", pos=1) %in% names(brw.attr))) warning(paste0(sum(!extract_string(v.nm, "\\|", pos=1) %in% names(brw.attr)), " nodes didn't map to 'brw.attr'. Unmapped nodes will be given seed values of 0."))
      ind = match(names(brw.attr), extract_string(v.nm, "\\|", pos=1))
      vertex_attr(graph, brw.attr.nm, ind[!is.na(ind)]) = brw.attr[!is.na(ind)]
    }else if(!(is.character(brw.attr) && length(brw.attr) > 1)) stop("Incorrect input for 'brw.attr'.")

    # Check FUN argument
    if(is.null(FUN) || (is.character(FUN) && length(FUN) == 1) || (is.function(FUN) && length(FUN) == 1)){
      # Assign actual function according to names...
      if(!is.function(FUN)) FUN = get_default_function(FUN)
    }else stop("Incorrect input type for 'FUN'.")
    # Check FUN.params argument
    if(!(is.null(FUN.params) || (is.list(FUN.params) && !is.null(names(FUN.params))))) stop("Incorrect input type for 'FUN.params'.")
    # Apply FUN with FUN.params. Check function results
    # Only pass function args if non-NULL
    if(is.null(FUN.params)){
      seeds = do.call(FUN, list(igraph::vertex_attr(graph, data.name)))
    }else seeds = do.call(FUN, c(list(igraph::vertex_attr(graph, data.name)), FUN.params))
    seeds[is.na(seeds)] = 0
    Z = (seeds - mean(seeds)) / sd(seeds)
    if(!all.good(seeds)) stop("Applying 'FUN' to 'data' resulted in either negative or all zero values.")
    # Set seeds and Z as vattrs
    igraph::vertex_attr(graph, "seeds") = seeds
    igraph::vertex_attr(graph, "Z") = Z
  }
  adj_mat = igraph::as_adjacency_matrix(graph = graph, attr = "weight", sparse = TRUE)

  # If multiplex, there are further pre-processing steps to link layers within multiplex components and to link bipartite edges to these layers
  if(multiplex){
    # For create_multiplex, I need to get every layer as a separate component in a list of adj mats
    multi.layer = unique(V(graph)$node_type[grep("_", V(graph)$node_type)])
    multi.adj = vector("list", length(multi.layer)); names(multi.adj) = multi.layer
    for(i in seq_along(multi.adj)){
      id = which(get.type(rownames(adj_mat), 2) %in% multi.layer[i] & n.char(get.type(rownames(adj_mat), 2), "_") == 1)
      multi.adj[[i]] = adj_mat[id,id]
    }

    mc = create_multiplex(adj_matrix = multi.adj, multiplex.comps = multiplex.comps)

    if(heterogeneous){
      # Get pairs of node types (no layer info) with bipartite edges (e.g., protein-metabolite, protein-gene)
      # Determining bipartite nodes, i.e., nodes that only have edges with another node type. These nodes will not have layer info, even if its node type is multiplex.
      uniq.types = unique(get.type(V(graph)$name, 3))
      node.types.l = get.type(V(graph)$name, 2)
      node.pairs = c()
      bp = logical(ncol(adj_mat))
      for(i in 1:length(bp)){
        id.tmp = adj_mat@i[(adj_mat@p[i] + 1):adj_mat@p[i+1]] + 1 # node ids of non-zero elements in column i
        bp[i] = any(uniq.types[uniq.types != node.types.l[i]] %in% node.types.l[id.tmp])
        if(bp[i]) node.pairs = c(node.pairs, sort(c(node.types.l[i], uniq.types[uniq.types != node.types.l[i]])))
      }
      node.pairs = unique(matrix(node.pairs, byrow = T, ncol = 2))

      # For each pair, get each component as adj mat and get the bp edges as a rectangular matrix
      # Input into 'transform_bipartite' (Treat monoplex as multiplex with 1 layer) and Merge the adjacency matrices
      new.adj.mat = vector("list", nrow(node.pairs))
      for(i in seq_along(new.adj.mat)){
        adj.tmp = vector("list", 2); names(adj.tmp) = node.pairs[i,]
        for(j in 1:2){
          # Handling potential multiplex components
          if(node.pairs[i,j] %in% names(mc)){
            adj.tmp[[j]] = mc[[node.pairs[i,j]]]$adj
          }else adj.tmp[[j]] = adj_mat[get.type(rownames(adj_mat),2) == node.pairs[i,j], get.type(rownames(adj_mat),2) == node.pairs[i,j]]
        }

        # Linking between node types (by layer if appropriate) using bipartite edges
        nl.in = n_layers(adj.tmp[[1]])
        nl.out = n_layers(adj.tmp[[2]])
        idc = which(bp & node.types.l == node.pairs[i,node.pairs[i,] == names(adj.tmp)[1]])
        idr = which(bp & node.types.l == node.pairs[i, node.pairs[i,] == names(adj.tmp)[2]])
        bp.mat = adj_mat[idr,idc]
        for(ll in seq_len(nl.out)){
          for(kk in seq_len(nl.in)){
            kk.tmp = transform_bipartite(col.mat = get_layer(adj.tmp[[1]], kk), row.mat = get_layer(adj.tmp[[2]], ll), bp.mat = bp.mat)
            if(kk == 1) kk.mat = kk.tmp else kk.mat = cbind(kk.mat, kk.tmp)
          }
          if(ll == 1) mat.tmp = kk.mat else mat.tmp = rbind(mat.tmp, kk.mat)
        }
        new.adj.mat[[i]] = cbind(rbind(adj.tmp[[1]], mat.tmp), rbind(Matrix::t(mat.tmp), adj.tmp[[2]]))
      }

      # Aggregate into an edge list and melt into a single graph
      edge_list = lapply(new.adj.mat, function(x){
        x = igraph::graph_from_adjacency_matrix(adjmatrix = x, mode = "undirected", weighted = TRUE)
        el = igraph::as_edgelist(graph = x, names = TRUE)
        if("weight" %in% igraph::edge_attr_names(x)){
          el =  cbind(el, igraph::E(x)$weight)
        }else el = cbind(el, 1)
        el
      })
      for(i in seq_along(edge_list)){
        if(i == 1){
          el = edge_list[[i]]
        }else{
          el = rbind(el, edge_list[[i]])
        }
      }
      g.tmp = igraph::graph_from_edgelist(el = el[,1:2], directed = FALSE)
      igraph::E(g.tmp)$weight = as.numeric(el[,3])
      g.tmp = igraph::simplify(graph = g.tmp, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "median")

      # Map v.attrs to new graph
      v.attrs = igraph::vertex_attr_names(graph)
      for(i in seq_along(v.attrs)){
        igraph::vertex_attr(g.tmp, v.attrs[i]) = igraph::vertex_attr(graph, v.attrs[i])[match(V(g.tmp)$name, V(graph)$name)]
      }

      graph = g.tmp

      # rm(g.tmp, edge_list, el, new.adj.mat, adj.tmp, bp.mat, mat.tmp)
    }else{ # If multiplex & homogeneous
      # Map v.attrs to new graph
      v.attrs = igraph::vertex_attr_names(graph)
      for(i in seq_along(v.attrs)){
        igraph::vertex_attr(mc[[1]]$g, v.attrs[i]) = igraph::vertex_attr(graph, v.attrs[i])[match(V(mc[[1]]$g)$name, V(graph)$name)]
      }

      graph = mc[[1]]$g

      # rm(mc, multi.adj)
    }
  }
  if(!igraph::is_connected(graph)){
    if(lcc){
      graph = largest_connected_component(graph)
    }else warning("Integrated graph is disconnected.")
  }
  return(graph)
}

#' @title Get largest connected component of a graph
#'
#' @description
#' Get largest connected component of a graph
#'
#' @param g igraph object
#'
#' @returns igraph object of largest connected component.
#'
largest_connected_component = function(g){
  if(!igraph::is_connected(g)){
    comps = igraph::components(g)
    largest_comp_id = which.max(comps$csize)
    g = igraph::induced_subgraph(g, which(comps$membership == largest_comp_id))
  }
  g
}

#' @title Link common nodes between layers of multiplex.
#'
#' @description
#' Given a list of inter-layer adjacency matrices (possibly from multiple multiplex components), this creates edges between common nodes across layers and places inter- and intra-layer matrices of a component in the same matrix.
#'
#' @param adj_matrix list of adjacency matrices.
#' @param multiplex.comps character vector of names of unique components that are multiplex.
#' @param data.type character vector of names of vertex attributes to preserve
#' @param brw.attr.nm character scalar  of vertex attribute name for biased random walk.
#'
#' @returns A list of fused multiplex components
#'
create_multiplex = function(adj_matrix, multiplex.comps, data.type, brw.attr.nm){
  type.tmp = extract_string(names(adj_matrix), "_", pos=1)
  multiplex.combined = vector("list", length(multiplex.comps)); names(multiplex.combined) = multiplex.comps
  for(i in seq_along(multiplex.combined)){
    # Identify IDs of layers
    l.ids = which(type.tmp %in% multiplex.comps[i])
    L = length(l.ids)
    layer.adjmat = adj_matrix[l.ids]
    new.mat = layer.adjmat[[1]]
    for(j in 2:L){
      for(ii in 1:(j-1)){
        # Construct inter-layer adj mat between layer j and layer ii
        col.nm = extract_string(colnames(layer.adjmat[[ii]])[Matrix::colSums(layer.adjmat[[ii]]) != 0], "\\|", pos=1)
        row.nm = extract_string(rownames(layer.adjmat[[j]])[Matrix::rowSums(layer.adjmat[[j]]) != 0], "\\|", pos=1)
        common = intersect(col.nm, row.nm)
        col.ids = match(paste(common, names(layer.adjmat)[ii], sep = "|"), colnames(layer.adjmat[[ii]]))
        row.ids = match(paste(common, names(layer.adjmat)[j], sep = "|"), rownames(layer.adjmat[[j]]))
        row.col.ids = matrix(c(row.ids, col.ids), ncol = 2)
        diag.tmp = Matrix::Matrix(0, nrow = nrow(layer.adjmat[[j]]), ncol = ncol(layer.adjmat[[ii]]), sparse = TRUE)
        diag.tmp[row.col.ids] = 1
        dimnames(diag.tmp) = list(rownames(layer.adjmat[[j]]), colnames(layer.adjmat[[ii]]))
        if(ii == 1){
          mat.tmp = diag.tmp
        }else{
          mat.tmp = cbind(mat.tmp, diag.tmp)
        }
      }
      new.mat = cbind(rbind(new.mat, mat.tmp), rbind(Matrix::t(mat.tmp), layer.adjmat[[j]]))
    }

    # Create multiplex graph. Will be disconnected if one layer has no nodes in common with another.
    # This is a problem in multiplex-homogeneous context, but not necessarily in a multiplex-heterogeneous context.
    multiplex.graph = igraph::graph_from_adjacency_matrix(adjmatrix = new.mat, mode = "undirected", weighted = TRUE)
    if(!igraph::is_connected(multiplex.graph)) stop(paste0("A layer in graph type \'", multiplex.comps[i],"\' has no nodes in common with any other layer."))

    # Add node_type to multiplex graph
    V(multiplex.graph)$node_type = extract_string(V(multiplex.graph)$name, "\\|", pos=2)

    multiplex.combined[[i]] = list(g = multiplex.graph, adj = new.mat)
  } # End For loop over 'multiplex.comps'
  # NB: Node types will be in blocks...

  return(multiplex.combined)
}

#' @title Transform a bipartite graph to conform column-wise to a matrix and row-wise to another matrix.
#'
#' @description
#' Transform bipartite graph _bp.mat_ to conform column-wise to _col.mat_ and row-wise to _row.mat_
#'
#' @param col.mat matrix to which _bp.mat_ is matched column-wise
#' @param row.mat matrix to which _bp.mat_ is matched row-wise
#' @param bp.mat bipartite matrix
#'
#' @returns a matrix
#'
transform_bipartite = function(col.mat, row.mat, bp.mat){
  r.t = extract_string(rownames(row.mat)[1], "\\|", pos=2)
  c.t = extract_string(rownames(col.mat)[1], "\\|", pos=2)

  col.nm = extract_string(colnames(col.mat), "\\|", pos=1)
  row.nm = extract_string(rownames(row.mat), "\\|", pos=1)

  if(all(grepl("\\|", rownames(bp.mat)))) rownames(bp.mat) = get.type(rownames(bp.mat),1)
  bp.row.nm = rownames(bp.mat)

  if(all(grepl("\\|", colnames(bp.mat)))) colnames(bp.mat) = get.type(colnames(bp.mat),1)
  bp.col.nm = colnames(bp.mat)

  c.id = which(bp.col.nm %in% col.nm)
  r.id = which(bp.row.nm %in% row.nm)
  bp.mat = bp.mat[r.id, c.id]

  missing.col.nm = col.nm[!col.nm %in% bp.col.nm]
  missing.row.nm = row.nm[!row.nm %in% bp.row.nm]

  tmat1 = rbind(bp.mat, Matrix::Matrix(0, nrow = length(missing.row.nm), ncol = ncol(bp.mat), sparse = TRUE))
  tmat2 = rbind(Matrix::Matrix(0, nrow = nrow(bp.mat), ncol = length(missing.col.nm), sparse = TRUE), Matrix::Matrix(0, nrow = length(missing.row.nm), ncol = length(missing.col.nm), sparse = TRUE))
  tmat = cbind(tmat1, tmat2)
  dimnames(tmat) = list(paste(c(bp.row.nm[r.id], missing.row.nm), r.t, sep = "|"), paste(c(bp.col.nm[c.id], missing.col.nm), c.t, sep = "|"))

  col.ids = match(colnames(col.mat), colnames(tmat))
  row.ids = match(rownames(row.mat), rownames(tmat))
  tmat = tmat[row.ids, col.ids]
  return(tmat)
}

#' @title Get a layer of a mono/multiplex component.
#'
#' @description
#' Get _l_th layer (in order of appearance) of a mono/multiplex component
#'
#' @param x matrix of a mono/multiplex component
#' @param l index of the layer to retrieve
#'
#' @returns matrix
#'
get_layer = function(x, l){
  t1 = extract_string(rownames(x), "\\|", pos=2)
  t2 = extract_string(t1, "_", pos=2)
  if(any(is.na(t2))){ # If there is nothing after "_", implying that this is a monoplex
    id = 1:length(t2)
  }else id = which(t2 == unique(t2)[l]) # layer l is the lth to appear in adjacency matrix
  x[id,id]
}

#' @title Get the number of layers for a mono/multiplex component
#'
#' @description
#' Get the number of layers for a mono/multiplex component
#'
#' @param x matrix
#'
#' @returns integer
#'
n_layers = function(x){
  t1 = unlist(lapply(strsplit(rownames(x), "\\|"), function(y) y[2]))
  t2 = unlist(lapply(strsplit(t1, "_"), function(y) y[2]))
  if(any(is.na(t2))) 1 else length(unique(t2))
}

#' @title Get a default function
#'
#' @description
#' Return a default function used to compute seed values for random walk with restart.
#'
#' For the input _f_, NULL means no transformation is done. 'binary' assumes that data is 0 or 1 so no transformation is done. 'shift_scale' is for data types with a range centered about zero and takes two parameters: _DOI_ (direction of interest: 1 for positive, -1 for negative, and 0 for both) and _w_ (numeric value between 0 and 1). It takes the absolute values of the data and then down-weights nodes that weren't in DOI by _w_ (if _DOI_=0, _w_ is coerced to 1). 'p_value' assumes the data are p-values and calculates _-log10()_ . 'exp' assumes the data are log fold changes and has a _DOI_ arg. For _DOI_=-1,1, it exponentiates the product of _DOI_ and _logFC_. For _DOI_=0, it exponentiates abs( _logFC_).
#'
#' @param f character string: 'binary', 'shift_scale', 'p_value', 'exp'. or NULL
#'
#' @returns a function
#'
get_default_function = function(f){
  if(is.null(f) || f == "binary"){
    res = function(x, ...) x
  }else if(f == "shift_scale"){
    res = function(x, DOI = 1, w = 0.5, ...){
      if(DOI == -1){
        ifelse(x < 0, abs(x), w * x)
      }else if(DOI == 1){
        ifelse(x > 0, x, w * abs(x))
      }else if(DOI == 0){
        abs(x)
      }else stop("DOI must be -1, 0, or 1.")
    }
  }else if(f == "p_value"){
    res = function(x, ...) -log(x + 1e-6, base = 10)
  }else if(f == "exp"){
    res = function(x, DOI = 1, ...){
      if(DOI == 0){
        exp(abs(x))
      }else if(DOI %in% c(-1,1)){
        exp(x * DOI)
      }else stop("DOI must be -1, 0, or 1.")
    }
  }else stop(paste0("Unrecognized character string \'", f,"\' for \'FUN\'."))
  res
}

#' @title Check if values of a vector are non-negative with at least one positive value.
#'
#' @description
#' Check if values of a vector are non-negative with at least one positive value.
#'
#' @param x numeric vector
#'
#' @returns logical: if all values are non-negative with at least one positive value.
#'
all.good = function(x) all(x >= 0) && any(x > 0)

#' @title The number of a user-specified string in a character string.
#'
#' @description
#' This returns the number of instances of _ch_ found in character string or vector _x_.
#'
#' @param x character scalar or vector to search
#' @param ch character to search for
#'
#' @returns integer
#'
n.char = function(x, ch) unlist(lapply(strsplit(x, ch), length)) - 1

#' @title Check the list names and node labels of network list inputs
#'
#' @description
#' Check the list names and node labels of network list inputs to ensure that they follow proper naming conventions. The algorithm needs to be able to associate each node in an input network with a component and layer. The component and layer information must be present as tags in the node labels (following 'name|component_layer' format) if it cannot be unambiguously determined for each node based on the network input list name.
#'
#' For example, a list input with the name 'gene' can either contain a monoplex component with or without '|gene' appended to node names, or a multiplex component with '|gene_layer' (e.g., 'gene_1', 'gene_2') appended to node names. Additionally, a multiplex component can be represented by giving each layer separately in the network input list, with possible names 'gene_1', 'gene_2', etc. Bipartite list inputs (i.e., containing more than one component) must be named with components separated by ';' (e.g., 'gene;drug').
#'
#' @details
#' For internal use within `create_integrated_graph()` only.
#'
#'
#' @param nn A named list of node labels
#' @param multiplex Logical. Is the integrated graph to be multiplex?
#' @param heterogeneous Logical. Is the integrated graph to be heterogeneous?
#'
#' @returns named list of updated node labels
#'
#' @examples
#' # Assume a multiplex-heterogeneous network of genes and drugs. The gene component is multiplex.
#' # Create hypothetical list of node names.
#' # The list names are taken from the network list inputs, and the node labels are taken from
#' #   the nodes of the network objects themselves.
#' node_names1 = list(gene = c("g1|gene_1", "g2|gene_1", "g3|gene_1", "g4|gene_1",
#'                             "g1|gene_2", "g2|gene_2", "g3|gene_2"),
#'                    drug = c("d1", "d2", "d3"),
#'                    "gene;drug" = c("d1|drug", "g1", "g2"))
#'
#' node_names2 = list(gene_1 = c("g1", "g2", "g3", "g4"),
#'                    gene_2 = c("g1", "g2", "g3"),
#'                    drug = c("d1", "d2", "d3"),
#'                    "gene;drug" = c("d1|drug", "d2", "g1|gene_1", "g1|gene_2", "g2|gene_1"))
#'
#' update_names1 = list_input_checks(nn = node_names1, multiplex = TRUE, heterogeneous = TRUE)
#' update_names2 = list_input_checks(nn = node_names2, multiplex = TRUE, heterogeneous = TRUE)
#' update_names1
#' update_names2
#'
#' @export
list_input_checks = function(nn, multiplex, heterogeneous){
  # 'nn' should be a named list of character vectors.
  error.message = "Incorrect labeling for graph components and node names. See function documentation for correct naming conventions."

  if(!multiplex){
    for(i in seq_along(nn)){
      if(length(unique(nn[[i]])) < length(nn[[i]])) stop("Not all node names are unique for certain network inputs.")
    }
  }else{
    if(!heterogeneous){ # Multiplex and NOT heterogeneous

      for(i in seq_along(nn)){ # For each input graph element
        if(!grepl(";", names(nn)[i])){ # If element i is not bipartite
          if(n.char(names(nn)[i], "_") == 1){ # If element i is a layer
            if(all(n.char(nn[[i]], "\\|") == 1) && length(unique(extract_string(nn[[i]], "\\|", pos=2))) == 1 && unique(extract_string(nn[[i]], "\\|", pos=2)) == names(nn)[i]){
              # GOOD. All names should be appended with '|component_layer'. Node names should not contain '|' apart from the delimiter. Appendage should match component name.
            }else if(all(!grepl("\\|", nn[[i]]))){ # If the labels don't include '|component_layer', add it.
              nn[[i]] = paste(nn[[i]], names(nn)[i], sep = "|")
            }else stop(error.message)
          }else if(n.char(names(nn)[i], "_") == 0){ # If element i is not explicitly labeled as multiplex
            if(all(n.char(nn[[i]], "\\|") == 1)){ # All labels contain '|'
              if(!(all(n.char(nn[[i]], "_") == 1) && length(unique(extract_string(extract_string(nn[[i]], "\\|", pos=2), "_", pos=1))) == 1 && unique(extract_string(extract_string(nn[[i]], "\\|", pos=2), "_", pos=1)) == names(nn)[i])){
                stop(error.message)
              }
            }else stop(error.message)
          }else stop(error.message)
        }else{ # if ';' is present, representing separate layers of multiplex
          tmp.nm = unlist(strsplit(names(nn)[i], ";"))
          if(!(all(n.char(tmp.nm, "_") == 1) && length(unique(extract_string(tmp.nm, "_", pos=1))) == 1 && length(unique(tmp.nm)) == length(tmp.nm))) stop(error.message)
          if(!(all(n.char(nn[[i]], "\\|") == 1) && all(extract_string(nn[[i]], "\\|", pos=2) %in% tmp.nm))) stop(error.message)
        }
      }

    }else{ # Multiplex and heterogeneous

      for(i in seq_along(nn)){
        if(!grepl(";", names(nn)[i])){ # If its mono/multiplex, not bipartite

          if(n.char(names(nn)[i], "_") == 1){ ## If one '_' is present, representing multiplex layer
            if(all(n.char(nn[[i]], "\\|") == 1) && length(unique(extract_string(nn[[i]], "\\|", pos=2))) == 1 && unique(extract_string(nn[[i]], "\\|", pos=2)) == names(nn)[i]){
              # GOOD. If they all have one '|', with one unique name after '|', and this unique name matches the name of the list.
            }else if(all(!grepl("\\|", nn[[i]]))){ # If no node name has a '|'
              nn[[i]] = paste(nn[[i]], names(nn)[i], sep = "|")
            }else stop(error.message)
          }else if(n.char(names(nn)[i], "_") == 0){ ## Else if NO '_' are present
            if(all(n.char(nn[[i]], "\\|") == 1)){
              if((all(n.char(nn[[i]], "_") == 1) || all(n.char(nn[[i]], "_") == 0)) && length(unique(extract_string(extract_string(nn[[i]], "\\|", pos=2), "_", pos=1))) == 1 && unique(extract_string(extract_string(nn[[i]], "\\|", pos=2), "_", pos=1)) == names(nn)[i]){
                # GOOD
              }else stop(error.message)
            }else if(all(!grepl("\\|", nn[[i]]))){ # If there is no "|" in any node name
              nn[[i]] = paste(nn[[i]], names(nn)[i], sep = "|") # Append the node type as given by list name.
            }else stop(error.message)
          }else stop(error.message)

        }else{ # If ";" is present, assume that this is a bipartite, heterogeneous graph and/or multiplex (e.g., mono vs mono, multi vs mono, multi vs multi, just multi)

          tmp.nm = unlist(strsplit(names(nn)[i], ";")) # Component names in bipartite graph
          if(any(grepl("_", tmp.nm))){ # If any designate a multiplex layer
            if(all(n.char(tmp.nm, "_") %in% c(1,0))){ # There can be at most one '_' per 'component_layer' term in list name
              multi.nm = tmp.nm[grep("_", tmp.nm)] # Names that contain '_'
              if(!all(n.char(nn[[i]], "\\|") %in% c(1,0))) stop(error.message) # There can be at most one '|' per node name, following the 'name' or 'name|component_layer' labeling format.
              if(!all(extract_string(nn[[i]][grep("\\|", nn[[i]])], "\\|", pos=2) %in% tmp.nm)) stop(error.message) # All components given in 'component_layer' must be in bipartite name of graph list input (with different components separated by ';')
              if(!all(tmp.nm %in% extract_string(nn[[i]][grep("\\|", nn[[i]])], "\\|", pos=2))) warning("Some node types designated in list names aren't present.") # All component_layers given in list name must be present in node labels, following 'name|component_layer' format.
            }else stop(error.message)
          }else{
            if(!(all(n.char(nn[[i]], "\\|") %in% c(1,0)))) stop(error.message) # There can be at most one '|' per node name, following the 'name' or 'name|component_layer' labeling format.
            id = grep("\\|", nn[[i]]) # ids of node labels with '|'
            tl = extract_string(nn[[i]][id], "\\|", pos=2) # Get the ''component_layer' info from these node labels
            # IF these are not all in list name AND not all component (not layer) info from node labels is present in list name, THEN error
            if(!all(tl %in% tmp.nm) && !all(extract_string(tl, "_", pos=1) %in% tmp.nm)) stop(error.message)
            # Scenario 1: all(component_layer %in% tmp.nm) = TRUE, all(component %in% tmp.nm) = FALSE
            #   If a node has component_layer info, and that is reflected in list name, that is sufficient. This is OK.
            # Scenario 2: all(component_layer %in% tmp.nm) = FALSE, all(component %in% tmp.nm) = TRUE
            #   If the component info is reflected in list name, but the component_layer info is not, that is OK.
            # Error wouldn't trip in both scenarios.
            # Nodes that don't have component_layer tags and don't map to anything given by node_type arg in create_integrated_graph() will be removed later on in that function.
          }

        }
      }
    }
  }
  if(!all(unlist(lapply(nn, function(x) length(unique(x)) == length(x))))) stop(error.message)
  nn
}

#' @title Perform checks on adjacency matrix inputs
#'
#' @description
#' Perform checks on adjacency matrix input to ensure proper data format.
#'
#' @param x adjacency matrix
#'
#' @returns adjacency matrix
#'
adj_mat_checks = function(x){
  if(!any(c("matrix", "array", "dgCMatrix", "dgRMatrix", "dgeMatirx", "dsCMatrix", "dsRMatrix", "dsyMatirx") %in% class(x))) stop("Unrecognized input type for adj_matrix.")
  if(is.null(dimnames(x))) stop("dimnames of adj_matrix is NULL.")
  # Convert to dgCMatrix
  x = methods::as(methods::as(methods::as(x, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  # If not a square matrix, make it square
  if(nrow(x) != ncol(x)){
    if(is.null(rownames(x)) || is.null(colnames(x))) stop("A non-square adjacency matrix is missing either row or column names.")
    rn.missing = setdiff(colnames(x), rownames(x))
    cn.missing = setdiff(rownames(x), colnames(x))
    nrn = c(rownames(x), rn.missing)
    ncn = c(colnames(x), cn.missing)
    x = cbind(rbind(x, Matrix::Matrix(0, nrow = length(rn.missing), ncol = ncol(x), sparse = T)),
               rbind(Matrix::Matrix(0, nrow = nrow(x), ncol = length(cn.missing), sparse = T), Matrix::Matrix(0, nrow = length(rn.missing), ncol = length(cn.missing), sparse = T)))
    dimnames(x) = list(nrn, ncn)
    x = x[order(nrn), order(ncn)]
  }
  if(is.null(rownames(x)) && !is.null(colnames(x))) rownames(x) = colnames(x)
  if(!is.null(rownames(x)) && is.null(colnames(x))) colnames(x) = rownames(x)
  if(any(!rownames(x) %in% colnames(x)) || any(!colnames(x) %in% rownames(x))) stop("Row and column names differ. If input is a square matrix, row and column names must match.")
  x[order(rownames(x)), order(colnames(x))]
}

#' @title Perform checks on edge list inputs
#'
#' @description
#' Perform checks on edge list input to ensure proper data format.
#'
#' @param x edge list in form of matrix, data.frame, or data.table
#'
#' @returns edge list
#'
el_checks = function(x){
  if(!any(c("matrix", "array", "dgeMatirx", "data.frame", "data.table") %in% class(x))) stop("Unrecognized input type for edge_list.")
  x = as.matrix(x)
  if(any(is.na(x))){
    warning("Removing NA values in Edge list.")
    x = stats::na.omit(x)
  }
  if(ncol(x) == 2){
    x = cbind(x, 1)
  }else if(ncol(x) < 2) stop("Edge lists must have at least two columns.")
  x
}

#' @title Perform checks on igraph inputs
#'
#' @description
#' Perform checks on igraph input to ensure proper data format.
#'
#' @param x igraph
#'
#' @returns igraph
#'
graph_checks = function(x){
  if(!"name" %in% igraph::vertex_attr_names(x)) stop("No 'name' vertex attribute.")
  nm = V(x)$name
  if("node_type" %in% igraph::vertex_attr_names(x)){
    id = grep("\\|", nm)
    if(length(id) != 0){
      if(any(n.char(nm[id], "\\|") != 1) || any(extract_string(nm[id], "\\|", pos=2) != V(x)$node_type[id])){
        stop("Incorrect node labels. See function documentation for correct naming conventions.")
      } #else warning("")
    }
    V(x)$name = paste(extract_string(nm, "\\|", pos=1), V(x)$node_type, sep = "|")
  }
  if(length(unique(V(x)$name)) < igraph::vcount(x)) stop("Not all node name/type combinations are unique.")
  x
}


#' @title Melt a graph
#'
#' @description
#' Melting a graph is finding nodes with the same node name and fusing them together so that there is just one node that is adjacent to all of the edges.
#'
#' @details
#' Fuses any nodes with same node name together, conserving any vertex attributes that may be present.
#'
#' @param g igraph
#' @param agg.method A function or NULL. The aggregation method for vertex attributes. Default NULL takes the first value.
#'
#' @returns igraph with vertex attributes preserved.
#'
melt_graph = function(g, agg.method = NULL){
  ## Collect any vertex attributes from graphs
  v.attrs = igraph::vertex_attr_names(g)
  if(length(v.attrs) != 0){
    tmp.scores = vector("list", length(v.attrs)); names(tmp.scores) = v.attrs
    for(i in seq_along(v.attrs)){
      tmp.scores[[i]] = igraph::vertex_attr(g, v.attrs[i])
      names(tmp.scores[[i]]) = V(g)$name
    }
  }
  edge_list = igraph::as_edgelist(graph = g, names = TRUE)
  if("weight" %in% igraph::edge_attr_names(g)){
    edge_list =  cbind(edge_list, igraph::E(g)$weight)
  }else edge_list = cbind(edge_list, 1)
  graph = igraph::graph_from_edgelist(el = edge_list[,1:2], directed = FALSE)
  igraph::E(graph)$weight = as.numeric(edge_list[,3])
  graph = igraph::simplify(graph = graph, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "median")
  if(length(v.attrs) != 0){
    # Redistribute the vertex attributes to the integrated graph. Takes the first vertex attribute value for nodes with multiple values.
    if(is.null(agg.method)){
      for(i in seq_along(v.attrs)){
        igraph::vertex_attr(graph, v.attrs[i]) = unname(tmp.scores[[i]][match(V(graph)$name, names(tmp.scores[[i]]))])
        # Switch NAs to zero
        # igraph::vertex_attr(graph, v.attrs[i], which(is.na(igraph::vertex_attr(graph, v.attrs[i])))) = 0
      }
    }else{
      # Get aggregation function
      agg.fun = get_aggregate_method(x = agg.method)
      for(i in seq_along(v.attrs)){
        if(is.numeric(tmp.scores[[i]])){
          if(any(tmp.scores[[i]] < 0) && agg.method %in% c("gmean", "hmean")) agg.fun.tmp = mean else agg.fun.tmp = agg.fun
          tmp = stats::aggregate(x = tmp.scores[[i]], by = list(names(tmp.scores[[i]])), FUN = agg.fun.tmp)
          igraph::vertex_attr(graph, v.attrs[i]) = tmp$x[match(V(graph)$name, tmp$Group.1)]
        }else{ # If non-numeric, take the first value
          igraph::vertex_attr(graph, v.attrs[i]) = unname(tmp.scores[[i]][match(V(graph)$name, names(tmp.scores[[i]]))])
        }
      }
    }
  }
  graph
}
