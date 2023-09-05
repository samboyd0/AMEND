#' @importFrom igraph vcount V V<- E E<- vertex_attr vertex_attr<-
#' @importFrom stats sd quantile

#' @title Identify active modules from an interaction network
#'
#' @description Identifies an active module, i.e., subnetwork with large experimental values, thorugh an iterative optimization procedure.
#'
#' @details
#'
#' Given an interaction network and experimental scores, AMEND attempts to find an active module by iteratively applying random walk with restart (RWR) and a heurstic solution to a maximum-weight connected subgraph problem.
#'
#' Briefly, RWR is applied to the network to get propagation scores which are shifted by some percentile (filtering rate, which decreases exponentially each iteration). The network with propagation scores are input into `heinz()`, an algorithm that heuristically finds a subnetwork with maximum node weight. This subnetwork is assigned a score based on its mean experimental values and mean core-clustering coefficients and then input into the next iteration. The algorithm stops when there is no change between iterations, and the subnetwork with maximum score is returned.
#'
#' This function can also accommodate heterogeneous networks with two node types. This information is incorporated in RWR for heterogeneous graphs (RWRH) proposed by [Li et al.](https://doi.org/10.1093/bioinformatics/btq108). This introduces the 'jump' and 'net1.weight' parameters.
#'
#' Since RWR requires seed values between 0 and 1, node scores often need to be transformed. 'seed.scheme' is used when data.type=ECI. 'zero_bottom' gives ECI values of zero the lowest weight, regardless of DOI. 'zero_middle' linearly transforms ECI values s.t. values in DOI are largest, values _not_ in DOI are smallest, and ECIs of zero are in the middle.
#'
#' The 'degree.bias' argument indicates whether degree bias in the network, which arises from technical and study biases, should be mitigated. This is done by scaling the transition matrix to be approximately bistochastic.
#'
#' See \url{link/to/paper} for more details on the AMEND algorithm. There have been slight modifications to the algorithm as presented in the original paper.
#'
#' @param graph igraph object with experimental data as vertex attribute. Either graph or adj_matrix and node_scores must be specified
#' @param adj_matrix Adjacency matrix of network. Either graph or adj_matrix and node_scores must be specified. Must have either colnames or rownames
#' @param node_scores Named vector of node scores
#' @param node_type Named vector of node types. Only considered if heterogeneous=TRUE. Can only accommodate two unique node types.
#' @param n Approximate size of the final module
#' @param data.type Type of input data. If graph is specified, must have matching vertex attribute
#' @param DOI Direction of Interest for experimental values. One of negative, positive, or both. When "both" is specified, the absolute value is taken. Only relevant for data type ECI and logFC
#' @param heterogeneous Logical. If TRUE, network is considered heterogeneous (two distinct node types, e.g., proteins and metabolites). If TRUE, node_type must be included as an argument or graph vertex attribute
#' @param normalize Normalization scheme of adjacency matrix for random walk with restart
#' @param eta Starting filtering rate. If NULL (default), a value is chosen based on the input network size and parameter 'n'
#' @param seed.weight Relative weight to give to nodes not in the DOI in random walk with restart, between 0 and 1
#' @param seed.scheme Determines how the feature-wise experimental values are transformed for use as seed values in RWR. See Details.
#' @param jump Probability of a random walker jumping from one network to the other in RWR (when heterogeneous=TRUE)
#' @param net1.weight Relative weight to give to nodes of network 1 (the first type to appear in node_type).
#' @param verbose Logical. Whether to output current iteration number to show progress
#' @param degree.bias Logical. Whether to mitigate degree bias through bistochastic scaling of transition matrix prior to RWR. See Details
#' @param identifier For use when performing many runs of AMEND to keep track of progress
#'
#' @return a named list with the following elements:
#' * module: the final module (i.e., subnetwork)
#' * score: final module score
#' * subnetworks: a list of node names contained in intermediate subnetworks
#' * stats: network statistics
#' * time: run time
#' * input_params: list of input parameters
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
#' # Use run_AMEND() with an igraph object with a vertex attribute matching data.type arg
#' subnet1 = run_AMEND(graph = glut4_graph, data.type = "ECI",
#'                 DOI = "negative")
#'
#' # Use run_AMEND() with an adjacency matrix and a vector of node scores
#' subnet2 = run_AMEND(adj_matrix = glut4_adjM, node_scores = eci_scores,
#'                 data.type = "ECI", DOI = "negative")
#' }
#'
#' @export
run_AMEND <- function(graph = NULL, adj_matrix = NULL, node_scores = NULL, node_type = NULL, n = 25,
                      data.type = c("ECI", "logFC", "p_val", "binary", "other"), DOI = c("positive", "negative", "both"), heterogeneous = FALSE,
                      normalize = c("degree", "core"), eta = NULL,
                      seed.weight = 0.5, seed.scheme = c("zero_bottom", "zero_middle"), jump = 0.5, net1.weight = 0.5,
                      verbose = TRUE, degree.bias = FALSE, identifier = 1){
  start_time = Sys.time()

  # Function for setting the starting filtering rate
  get.eta0 = function(g, n, method = c("difference", "ratio", "log.ratio")){
    # eta is linearly increasing with D(N,n), where D is some distance metric between the input network size N and the final subnetwork size n
    # D(N,n) = N - n; for method = "difference"
    # D(N,n) = N / n; for method = "ratio"
    # D(N,n) = log(N / n); for method = "log.ratio"
    # Boundaries within which eta changes linearly with D(N,n)
    # N: [500, 15000]
    # n: [5, 200]
    # eta: [0.3, 0.8]
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

  # Variables not set by user, but could change in future
  ma_window = 5 # Moving average window... set to Inf for a running average

  data.type = match.arg(data.type)
  DOI = match.arg(DOI)
  normalize = match.arg(normalize)
  seed.scheme = match.arg(seed.scheme)

  # Ensuring adj_matrix is the correct class
  if(all(class(adj_matrix) != "matrix") && !is.null(adj_matrix)){
    adj_matrix = as.matrix(adj_matrix)
  }

  # NB: If igraph object is given, we want to preserve any node/edge attrs it may have, so igraph takes priority over adj_matrix
  if(!is.null(graph)){
    # Check name attr
    if(!"name" %in% igraph::vertex_attr_names(graph)) stop("No \'name\' vertex attribute.")
    # Check node scores
    if(!data.type %in% igraph::vertex_attr_names(graph)){
      if(is.null(node_scores)){
        stop("There is no vertex attribute \'", data.type, "\' or named vector of node scores. At least one must be given")
      }else{
        if(is.null(names(node_scores))) stop("node_scores is not a named vector")
        if(any(!V(graph)$name %in% names(node_scores))) stop("node_scores does not contain all nodes in graph")
        ind = match(names(node_scores), V(graph)$name)
        vertex_attr(graph, data.type, ind[!is.na(ind)]) = node_scores[!is.na(ind)]
      }
    }
    # Check node type
    if(heterogeneous){
      if(!"node_type" %in% igraph::vertex_attr_names(graph)){
        if(is.null(node_type)){
          stop("There is no vertex attribute \'node_type\' or named vector of node types. At least one must be given for heterogeneous graphs.")
        }else{
          if(is.null(names(node_type))) stop("node_type is not a named vector")
          if(any(!V(graph)$name %in% names(node_type))) stop("node_type does not contain all nodes in graph.")
          ind = match(names(node_type), V(graph)$name)
          vertex_attr(graph, "node_type", ind[!is.na(ind)]) = node_type[!is.na(ind)]
        }
      }
      if(length(unique(V(graph)$node_type)) > 2) stop("Number of unique node types is greater than 2. Currently, only heterogeneous graphs of two node types is supported.")
    }
    # Check edge weights... if igraph is given w/o edge attr, assume it is unweighted, even if a weighted adj matrix is given
    if(!"weight" %in% igraph::edge_attr_names(graph)) E(graph)$weight = 1
  }else if(is.null(adj_matrix) || is.null(node_scores)){
    if(heterogeneous && is.null(node_type)){
      stop("No inputs for graph, adj_matrix, node_scores, or node_type (for heterogeneous graphs).")
    }else stop("No inputs for graph, adj_matrix, or node_scores.")
  }else{
    if(is.null(names(node_scores))) stop("node_scores is not a named vector.")
    if(is.null(rownames(adj_matrix)) && is.null(colnames(adj_matrix))){
      stop("adj_matrix does not have rownames or colnames. One of these must be present.")
    }else{
      if(is.null(rownames(adj_matrix))) rownames(adj_matrix) = colnames(adj_matrix)
      if(is.null(colnames(adj_matrix))) colnames(adj_matrix) = rownames(adj_matrix)
    }
    # Matching node_scores to order of rows in adj_matrix
    if(any(!rownames(adj_matrix) %in% names(node_scores))) stop("node_scores does not contain all nodes in adj_matrix.")
    ind = match(names(node_scores), rownames(adj_matrix))
    node_scores = node_scores[!is.na(ind)] # removes elements that didn't map
    node_scores[ind[!is.na(ind)]] = node_scores # Get correct order
    # Creating igraph object
    graph = igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = "weight", diag = FALSE)
    graph = igraph::set_vertex_attr(graph = graph, name = data.type, value = node_scores)
    # Matching node_type to order of rows in adj_matrix
    if(heterogeneous){
      if(is.null(node_type)){
        stop("node_type is not given.")
      }else if(is.null(names(node_type))) stop("node_type is not a named vector.")
      if(any(!rownames(adj_matrix) %in% names(node_type))) stop("node_type does not contain all nodes in adj_matrix.")
      ind = match(names(node_type), rownames(adj_matrix))
      node_type = node_type[!is.na(ind)] # removes elements that didn't map
      node_type[ind[!is.na(ind)]] = node_type # Get correct order
      graph = igraph::set_vertex_attr(graph = graph, name = "node_type", value = node_type) # set as vertex attr
      if(length(unique(node_type)) > 2) stop("Number of unique node types is greater than 2. Currently, only heterogeneous graphs of two node types is supported.")
    }
  }
  # At this point, evrything is in 'graph',
  # so the goal here is to get vertex sequence in graph s.t. all nodes of a given type are grouped together
  if(heterogeneous){
    nt = unique(V(graph)$node_type)
    p = Matrix::invPerm(c(which(V(graph)$node_type == nt[1]), which(V(graph)$node_type == nt[2])))
    graph = igraph::permute(graph, p)
  }
  # Transforming input data
  graph = data_preprocessing(graph, data.type, DOI, seed.weight, seed.scheme, heterogeneous)

  repeat{
    all_scores <- list()
    all_nets <- list()

    subg <- vector(mode = "list", length = 2)
    subg[[1]] <- graph
    k <- 1

    avg.rd = 0.07 # initial average rate difference

    if(is.null(eta)) eta = get.eta0(graph, n, "log.ratio")

    if(verbose) message(paste("Starting filtering rate:", round(eta, 4)))
    repeat{
      if(verbose) message(paste("Iteration:", k))
      if(k == 1){
        rate.difference = avg.rd
        e = eta
        decay.N = vcount(graph)

        adj_og <- igraph::as_adjacency_matrix(subg[[1]], type = "both", attr = "weight", sparse = T)
        adj <- adj_og
      }else{
        rate.difference = c(avg.rd, do.call(c, lapply(all_scores, function(x) x[10] - x[9])))
        rate.difference = rate.difference[ifelse(length(rate.difference) >= ma_window, length(rate.difference) - ma_window + 1, 1):length(rate.difference)]
        e = all_scores[[k-1]][9]

        if(k == 2){
          decay.N = vcount(graph)
        }else decay.N = length(all_nets[[k-2]])

        adj <- adj[-rm.id, -rm.id]
      }
      # Setting shifting percentile (i.e., filtering rate) for this iteration
      decay <- find_decay(N = decay.N, eta0 = e, n = n, rate.diff = rate.difference)
      l <- exp_filtering_rate(eta0 = e, d = decay)
      ll <- l[ifelse(k == 1, 1, 2)]

      # Determining whether graph is still heterogeneous, b/c all of one node type may have been filtered out
      if(heterogeneous){
        if(length(unique(V(subg[[1]])$node_type)) < 2){
          heterogeneous = FALSE
          if(verbose) message(paste0("All nodes of type \'", setdiff(nt, unique(V(subg[[1]])$node_type)), "\' have been filtered out."))
        }
      }
      # Normalize adjacency matrix to get transition matrix
      n.adjM <- transition_matrix(g = subg[[1]], adjM = adj, norm = normalize, heterogeneous = heterogeneous, node_type = V(subg[[1]])$node_type, jump = jump)

      # Perform bistochastic scaling to mitigate degree bias, if specified
      if(degree.bias) n.adjM = bistochastic_scaling(trans_mat = n.adjM)

      # Create matrix of seed values
      Seeds <- matrix(V(subg[[1]])$seeds, ncol = 1, dimnames = list(V(subg[[1]])$name))

      # Choosing Restart parameter value through a 'normal' grid search
      rgs <- restart_grid_search(subg[[1]], n.adjM, Seeds, ll, heterogeneous, V(subg[[1]])$node_type, net1.weight, k)
      subg[[2]] <- rgs[[1]]

      # break out of repeat loop if there is no change in network or if subnet size <= 2
      if(vcount(subg[[2]]) == vcount(subg[[1]]) || vcount(subg[[2]]) <= 2) break

      mCC <- rgs[[2]][1] # Mean core-clustering coefficient
      mZ <- rgs[[2]][2] # Mean standardized experimental scores
      n_nodes <- rgs[[2]][3] # Size of subnetwork

      sn_score <- mCC * mZ # subnetwork score

      all_scores[[length(all_scores) + 1]] <- c(decay, rgs[[3]], sn_score, mZ, mCC, n_nodes, igraph::ecount(subg[[2]]), igraph::edge_density(subg[[2]]), ll, mean(!V(subg[[1]])$name %in% V(subg[[2]])$name))
      all_nets[[length(all_nets) + 1]] = V(subg[[2]])$name

      # IDs of nodes that were removed this iteration
      rm.id <- which(!V(subg[[1]])$name %in% V(subg[[2]])$name)

      k <- k + 1
      subg[[1]] <- subg[[2]]
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
      if(verbose) message(paste0("Algorithm didn't converge with eta=", round(eta, 4), ". Trying a larger eta value. ID=", identifier))
      eta = eta * 1.1 # increase by 10%
    }
  }

  # Choosing subnetwork with maximum score and (in the case of ties) closest to final module size 'n'
  # and in case of two networks with same max score and same distance from 'n', choose larger subnetwork
  if(sum(all_scores[size.cond, 3] == max(all_scores[size.cond, 3])) > 1){ # Handling ties
    d = abs(all_scores[,6] - n)
    ids = which(all_scores[,3] == max(all_scores[size.cond, 3]) & size.cond)
    max.id = which(all_scores[,3] == max(all_scores[size.cond, 3]) & d == min(d[ids]) & size.cond)[1]
  }else max.id = which(all_scores[,3] == max(all_scores[size.cond, 3]) & size.cond)

  best_sn = igraph::induced_subgraph(graph, which(V(graph)$name %in% all_nets[[max.id]]))
  best_score = all_scores[max.id,3]

  end_time <- Sys.time()
  time <- end_time - start_time
  message(paste0("*** Converged! *** ID=", identifier))

  rm(list = ls()[!ls() %in% c("eta", "graph", "n", "data.type", "DOI", "adj_matrix", "normalize",
                              "seed.weight", "seed.scheme", "best_sn", "best_score", "all_scores", "time", "all_nets")])

  input_params = list(eta0 = eta, ig = graph, n = n, data.type = data.type, DOI = DOI, adjM = adj_matrix,
                      norm = normalize, seed.weight = seed.weight, seed.scheme = seed.scheme)

  return(list(module = best_sn, score = best_score, subnetworks = all_nets, stats = all_scores, time = time, input_params = input_params))
}

#' @title Pre-processing of experimental data for input into AMEND and RWR
#'
#' @description Pre-processing of experimental data for input into AMEND and RWR. For use inside `run_AMEND()`
#'
#' @details Transforms experimental data for use with Random walk with restart (RWR), which requires non-negative seed values with at least one non-zero element. For subnetwork scoring, the experimental scores need to be standardized w.r.t. all nodes in the input network.
#' For data.type "ECI", the seeding scheme involves taking the absolute value of the ECIs then weighting those nodes NOT in the Direction of Interest (DOI) by 'seed.weight'
#' For data.type "logFC", the log fold changes are multiplied by -1 if the DOI is negative, 1 if the DOI is positive, and the absolute value is taken if the DOI is both. Then the values are exponentiated
#' For data.type "p_val", the p-values are transformed by -log10(p-value)
#' For data.type "binary", weight is assigned uniformly to "active" nodes, i.e., nodes with a value of 1 in the binary vector.
#' For data.type "other", no transformation is performed.
#'
#' @param graph Input graph
#' @param data.type One of ECI, logFC, p_val, or binary. If graph is specified, must have matching vertex attribute
#' @param DOI Direction of Interest for experimental values. One of negative, positive, or both. When "both" is specified, the absolute value is taken. Only relevant for data type ECI and logFC
#' @param seed.weight Relative weight to give to nodes not in the DOI in random walk with restart, between 0 and 1
#' @param seed.scheme Determines how the feature-wise experimental values are transformed for use as seed values in RWR. See Details.
#' @param heterogeneous Logical. If TRUE, network is considered heterogeneous (two distinct node types, e.g., proteins and metabolites). If TRUE, node_type must be included as an argument or graph vertex attribute
#'
#' @return igraph object with added 'seeds' and 'Z' vertex attributes
#'
#' @examples
#' # Attach igraph package
#' library(igraph)
#'
#' # Inspect the igraph object included in AMEND package
#' # One can see it has vertex attributes name, symbol, and ECI
#' glut4_graph
#' plot(density(V(glut4_graph)$ECI))
#'
#' # Transform experimental values (ECI) for multiple values of 'seed.weight'
#' new_graph1 = AMEND:::data_preprocessing(graph = glut4_graph, data.type = "ECI",
#'                                         DOI = "negative", seed.weight = 0.5)
#' new_graph2 = AMEND:::data_preprocessing(graph = glut4_graph, data.type = "ECI",
#'                                         DOI = "negative", seed.weight = 0.2)
#'
#' # Comparing
#' plot(density(V(new_graph1)$seeds))
#' plot(density(V(new_graph2)$seeds))
#'
#' plot(V(new_graph1)$ECI, V(new_graph1)$Z, type = "b",
#'      xlim = c(-1, 1), xlab = "ECI", ylab = "Standardized ECI")
#'
data_preprocessing = function(graph, data.type = c("ECI", "logFC", "p_val", "binary", "other"), DOI = c("positive", "negative", "both"),
                              seed.weight = 0.5, seed.scheme = c("zero_bottom", "zero_middle"), heterogeneous = FALSE){
  data.type = match.arg(data.type)
  DOI = match.arg(DOI)
  seed.scheme = match.arg(seed.scheme)

  # Transforming input data
  if(data.type == "ECI"){
    if(seed.weight < 0){
      stop("Seed weight needs to be non-negative. Should be between 0 and 1.")
    }else if(seed.weight > 1){
      warning("Seed weight should be between 0 and 1.")
    }
    if(DOI == "both"){
      vertex_attr(graph, "seeds") = abs(V(graph)$ECI)
      vertex_attr(graph, "Z") = (V(graph)$seeds - mean(V(graph)$seeds)) / sd(V(graph)$seeds)
    }else{
      s <- ifelse(DOI == "positive", 1, -1)
      if(seed.scheme == "zero_bottom"){
        vertex_attr(graph, "seeds") = ifelse(s == rep(1, vcount(graph)), ifelse(V(graph)$ECI > 0, V(graph)$ECI, seed.weight * abs(V(graph)$ECI)), ifelse(V(graph)$ECI < 0, abs(V(graph)$ECI), seed.weight * V(graph)$ECI))
      }else if(seed.scheme == "zero_middle"){
        vertex_attr(graph, "seeds") <- s*V(graph)$ECI + 1
      }else stop("Unknown seed scheme specified.")
      vertex_attr(graph, "Z") <- (V(graph)$ECI - mean(V(graph)$ECI))/sd(V(graph)$ECI) * s
    }
  }else if(data.type == "logFC"){
    if(DOI == "both"){
      vertex_attr(graph, "seeds") <- exp(abs(V(graph)$logFC))
      vertex_attr(graph, "Z") <- (abs(V(graph)$logFC) - mean(abs(V(graph)$logFC)))/sd(abs(V(graph)$logFC))
    }else if(DOI %in% c("positive", "negative")){
      s <- ifelse(DOI == "positive", 1, -1)
      vertex_attr(graph, "seeds") <- exp(s*V(graph)$logFC)
      vertex_attr(graph, "Z") <- (V(graph)$logFC - mean(V(graph)$logFC))/sd(V(graph)$logFC) * s
    }
  }else if(data.type == "p_val"){
    vertex_attr(graph, "seeds") <- -log(vertex_attr(graph, "p_val") + 0.001, base = 10)
    vertex_attr(graph, "Z") <- (V(graph)$seeds - mean(V(graph)$seeds))/sd(V(graph)$seeds)
  }else if(data.type == "binary"){
    vertex_attr(graph, "seeds") = vertex_attr(graph, "binary") / sum(vertex_attr(graph, "binary"))
    vertex_attr(graph, "Z") = (V(graph)$seeds - mean(V(graph)$seeds))/sd(V(graph)$seeds)
  }else if(data.type == "other"){
    vertex_attr(graph, "seeds") = vertex_attr(graph, "other")
    vertex_attr(graph, "Z") = (V(graph)$seeds - mean(V(graph)$seeds))/sd(V(graph)$seeds)
  }else stop(paste0("Data type \'", data.type, "\' not recognized."))
  # Addressing possible NA values
  V(graph)$seeds[is.na(V(graph)$seeds)] = 0
  V(graph)$Z[is.na(V(graph)$Z)] = 0
  return(graph)
}

#' @title Grid search for the restart parameter in RWR
#'
#' @description Grid search for the restart parameter in RWR. For use inside `run_AMEND()`
#'
#' @details This is a simple grid search for the restart probability parameter in RWR.
#'
#' @param ig Input graph
#' @param n.adj.M Normalized adjacency matrix
#' @param seeds Vector of seed values
#' @param filtering_rate Quantile for shifting the raw RWR scores
#' @param heterogeneous Logical. If TRUE, network is considered heterogeneous (two distinct node types, e.g., proteins and metabolites). If TRUE, node_type must be included as an argument or graph vertex attribute
#' @param node_type Named vector of node types. Only considered if heterogeneous=TRUE. Can only accommodate two unique node types.
#' @param net1.weight Relative weight to give to nodes of network 1 (the first type to appear in node_type).
#' @param iteration Current iteraiton of the AMEND algorithm
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
#'
#' # Creating a named vector of seed values
#' seeds = runif(8)
#' names(seeds) = 1:8
#'
#' search = AMEND:::restart_grid_search(ig = g, n.adj.M = adj_norm,
#'                                      seeds = seeds, filtering_rate = 0.3)
#' plot(search[[1]])
#'
restart_grid_search <- function(ig, n.adj.M, seeds, filtering_rate, heterogeneous = FALSE, node_type = NULL, net1.weight = 0.5, iteration = 1){
  if(iteration == 1){
    grid <- seq(0.75, 0.95, by = 0.05)
  }else grid <- seq(0.5, 0.95, by = 0.05)

  score_metrics <- matrix(nrow = length(grid), ncol = 4)
  nets <- vector(mode = "list", length = length(grid))
  rm.id = c()
  for(i in seq_along(grid)){
    # RWR
    PTmatrix2 <- RWR(nadjM = n.adj.M, setSeeds = seeds, restart = grid[i], heterogeneous = heterogeneous, node_type = node_type, net1.weight = net1.weight)
    rwr_score <- as.vector(PTmatrix2)

    # Scores for Maximum Scoring Subgraph Algorithm
    rwr_score <- rwr_score - quantile(rwr_score, filtering_rate) # if filtering rate close to 0, quantile takes the minimum of rwr_score
    names(rwr_score) <- V(ig)$name
    nets[[i]] <- heinz(ig, rwr_score)

    n.v <- vcount(nets[[i]])
    mZ <- sum(V(nets[[i]])$Z) / n.v
    mCC = sum(core_cc(nets[[i]])) / n.v

    # To prevent filtering out too many nodes at one iteration
    # This is the difference between the observed and theoretical filtering rates
    if((vcount(ig) - n.v) / vcount(ig) - filtering_rate >= 0.25){ # 3/16/2023: used to be 0.5
      rm.id = c(rm.id, i)
      next
    }

    score_metrics[i,] <- c(mCC, mZ, n.v, grid[i])
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
    new.g = ig
    n.v = vcount(new.g)
    new.scores = c(sum(core_cc(new.g)) / n.v, sum(V(new.g)$Z) / n.v, n.v, NA)
    best.r = NA
  }

  return(list(new.g, new.scores, best.r))
}

#' @title Create a transition matrix from the adjacency matrix of a graph
#'
#' @description `transition_matrix()` creates a transition matrix from the adjacency matrix of a graph for input into random walk with restart (RWR).
#'
#' @details There are two normalization schemes available.
#' * degree: divides each column of the adjacency matrix by the column total
#' * core: uses the concept of node coreness. Details can be found in \url{https://academic.oup.com/nar/article/48/17/e98/5879427}
#'
#' @param g Input graph
#' @param adjM Adjacency matrix
#' @param norm Normalization method
#' @param heterogeneous Logical. If TRUE, network is considered heterogeneous (two distinct node types, e.g., proteins and metabolites). If TRUE, node_type must be included as an argument or graph vertex attribute.
#' @param node_type Named vector of node types. Only considered if heterogeneous=TRUE. Can only accommodate two unique node types.
#' @param jump Probability of a random walker jumping from one network to the other in RWR (when heterogeneous=TRUE).
#' @param k Value between 0 and 1 indicating by how much to penalize for node degree when calculating new adjacency matrix edge weights for the "modified_degree" option.
#'
#' @return transition matrix
#'
#' @examples
#'
#'
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
#'                0, 0, 0, 0, 1, 0, 0, 0), nrow = 8)
#' g = graph_from_adjacency_matrix(adjm, mode = 'undirected')
#'
#' adj_norm = AMEND:::transition_matrix(g, adjm, norm = "degree")
#'
#' @export
transition_matrix <- function(g, adjM, norm = c("degree", "modified_degree"), heterogeneous = FALSE, node_type = NULL, jump = 0.5, k = 0.5){

  #=== To-Do: Generalize heterogeneous code to an arbitrary number of different node types.

  sum2one <- function(X) {
    if(!"dgCMatrix" %in% class(X)) X = Matrix::Matrix(X, sparse = TRUE)
    inv_col_sum = Matrix::Diagonal(x = abs(Matrix::colSums(X))^(-1))
    res = X %*% inv_col_sum
    res[is.na(res)] = 0
    res
  }
  entropy = function(x){
    tmp = ifelse(round(x, 10) == 0, 0, -x * log(x))
    sum(tmp)
  }
  stationary.distr = function(x){
    e = Re(RSpectra::eigs(A = x, k = 1, which = "LM")$vectors[,1])
    tmp = e / sum(e)
    ifelse(tmp < 0, 0, tmp)
  }
  modify_adj_mat = function(adjM, k){
    # For symmetric matrices, inv.col.sum = inv.row.sum
    inv.col.sum <- Matrix::Diagonal(x = Matrix::colSums(adjM)^(-k))
    inv.row.sum <- Matrix::Diagonal(x = Matrix::rowSums(adjM)^(-k))
    adjM.mod <- inv.row.sum %*% adjM %*% inv.col.sum
    return(adjM.mod)
  }

  # Coerce to a general, column-compressed sparse matrix from Matrix package
  if(!"dgCMatrix" %in% class(adjM)) adjM = Matrix::Matrix(adjM, sparse = TRUE) # methods::as(adjM, "dgCMatrix")

  if(!heterogeneous){
    norm <- match.arg(norm)

    if(0){ # norm == "core"
      adj.dense = as.matrix(adjM)
      core <- sna::kcores(adj.dense, mode = "graph", ignore.eval = T)
      nadjM = matrix(0, nrow = nrow(adjM), ncol = nrow(adjM))
      for(j in 1:ncol(nadjM)){
        id.tmp = adjM@i[(adjM@p[j] + 1):adjM@p[j+1]] + 1 # row ids of non-zero elements in column j
        denom = sum(core[id.tmp])
        nadjM[id.tmp, j] = core[id.tmp] / denom
      }
      nadjM = methods::as(nadjM, "dgCMatrix")
    }else if(norm == "degree"){
      wt <- Matrix::Diagonal(x = Matrix::colSums(adjM)^(-1))
      nadjM <- adjM %*% wt
    }else if(norm == "modified_degree"){
      if(is.null(k)){
        node.strength = Matrix::colSums(adjM)
        k.grid = seq(0.5, 0.8, 0.05)
        ent = numeric(length(k.grid))
        for(i in seq_along(k.grid)){
          wt <- Matrix::Diagonal(x = Matrix::colSums(adjM)^(-k.grid[i]))
          adjM.mod <- wt %*% adjM %*% wt
          wt.mod <- Matrix::Diagonal(x = Matrix::colSums(adjM.mod)^(-1))
          nadjM <- adjM.mod %*% wt.mod
          sdist = stationary.distr(nadjM)
          ent[i] = entropy(sdist)
        }
        k.best = k.grid[which.max(ent)]
        wt <- Matrix::Diagonal(x = Matrix::colSums(adjM)^(-k.best))
        adjM.mod <- wt %*% adjM %*% wt
        wt.mod <- Matrix::Diagonal(x = Matrix::colSums(adjM.mod)^(-1))
        nadjM <- adjM.mod %*% wt.mod
      }else{
        wt <- Matrix::Diagonal(x = Matrix::colSums(adjM)^(-k))
        adjM.mod <- wt %*% adjM %*% wt
        wt.mod <- Matrix::Diagonal(x = Matrix::colSums(adjM.mod)^(-1))
        nadjM <- adjM.mod %*% wt.mod
      }
    }else{
      nadjM <- adjM
    }
    dimnames(nadjM) = dimnames(adjM)
    return(nadjM)
  }else{
    norm <- match.arg(norm)

    if(is.null(node_type)) stop("node_type is null.")

    # Unique node types
    nt = unique(node_type)
    if(length(nt) > 2) stop("There are more than 2 node types.")

    adj.1 = adjM[node_type == nt[1], node_type == nt[1]]
    adj.2 = adjM[node_type == nt[2], node_type == nt[2]]
    adj.12 = adjM[node_type == nt[1], node_type == nt[2]]
    adj.21 = Matrix::t(adj.12)

    # Checking row and column names
    if(is.null(rownames(adj.1))) rownames(adj.1) = rownames(adjM)[node_type == nt[1]]
    if(is.null(colnames(adj.1))) colnames(adj.1) = colnames(adjM)[node_type == nt[1]]

    if(is.null(rownames(adj.2))) rownames(adj.2) = rownames(adjM)[node_type == nt[2]]
    if(is.null(colnames(adj.2))) colnames(adj.2) = colnames(adjM)[node_type == nt[2]]

    if(is.null(rownames(adj.12))) rownames(adj.12) = rownames(adjM)[node_type == nt[1]]
    if(is.null(colnames(adj.12))) colnames(adj.12) = colnames(adjM)[node_type == nt[2]]

    if(is.null(rownames(adj.21))) rownames(adj.21) = rownames(adjM)[node_type == nt[2]]
    if(is.null(colnames(adj.21))) colnames(adj.21) = colnames(adjM)[node_type == nt[1]]

    # Determining bipartite nodes
    bp = logical(ncol(adjM))
    for(i in 1:length(bp)){
      id.tmp = adjM@i[(adjM@p[i] + 1):adjM@p[i+1]] + 1 # node ids of non-zero elements in column j
      bp[i] = all(nt %in% node_type[id.tmp])
    }
    bp.1 = bp[node_type == nt[1]]
    bp.2 = bp[node_type == nt[2]]

    # Modifying the adjacency matrices
    if(norm == "modified_degree"){
      if(is.null(k)) k = 0.5

      adj.1 = modify_adj_mat(adjM = adj.1, k = k)
      adj.2 = modify_adj_mat(adjM = adj.2, k = k)
      adj.12 = modify_adj_mat(adjM = adj.12, k = k)
      adj.21 = modify_adj_mat(adjM = adj.21, k = k)

      norm = "degree"
    }

    if(0){ # norm == "core"
      # core <- igraph::coreness(g)
      adj.dense = as.matrix(adjM)
      core <- sna::kcores(adj.dense, mode = "graph", ignore.eval = T)
      core.1 = core[node_type == nt[1]]
      core.2 = core[node_type == nt[2]]

      # PPIs
      nadj.1 = matrix(0, nrow = nrow(adj.1), ncol = nrow(adj.1), dimnames = dimnames(adj.1))
      c.sum = Matrix::colSums(adj.1)
      jump.tmp = ifelse(bp.1, 1 - jump, 1)
      for(j in which(c.sum != 0)){
        id.tmp = adj.1@i[(adj.1@p[j] + 1):adj.1@p[j+1]] + 1 # row ids of non-zero elements in column j
        denom = sum(core.1[id.tmp])
        nadj.1[id.tmp, j] = jump.tmp[j] * core.1[id.tmp] / denom
      }
      nadj.1 = methods::as(nadj.1, "dgCMatrix")

      # MMIs
      nadj.2 = matrix(0, nrow = nrow(adj.2), ncol = nrow(adj.2), dimnames = dimnames(adj.2))
      c.sum = Matrix::colSums(adj.2)
      jump.tmp = ifelse(bp.2, 1 - jump, 1)
      for(j in which(c.sum != 0)){
        id.tmp = adj.2@i[(adj.2@p[j] + 1):adj.2@p[j+1]] + 1 # row ids of non-zero elements in column j
        denom = sum(core.2[id.tmp])
        nadj.2[id.tmp, j] = jump.tmp[j] * core.2[id.tmp] / denom
      }
      nadj.2 = methods::as(nadj.2, "dgCMatrix")

      # pm: proteins in rows, metabolites in columns
      nadj.12 = matrix(0, nrow = nrow(adj.12), ncol = ncol(adj.12), dimnames = dimnames(adj.12))
      c.sum = Matrix::colSums(adj.12)
      for(j in which(c.sum != 0)){
        id.tmp = adj.12@i[(adj.12@p[j] + 1):adj.12@p[j+1]] + 1 # row ids of non-zero elements in column j
        denom = sum(core.1[id.tmp])
        nadj.12[id.tmp, j] <- jump * core.1[id.tmp] / denom
      }
      nadj.12 = methods::as(nadj.12, "dgCMatrix")

      # mp: metabolites in rows, proteins in columns
      nadj.21 = matrix(0, nrow = nrow(adj.21), ncol = ncol(adj.21), dimnames = dimnames(adj.21))
      c.sum = Matrix::colSums(adj.21)
      for(j in which(c.sum != 0)){
        id.tmp = adj.21@i[(adj.21@p[j] + 1):adj.21@p[j+1]] + 1 # row ids of non-zero elements in column j
        denom = sum(core.2[id.tmp])
        nadj.21[id.tmp, j] <- jump * core.2[id.tmp] / denom
      }
      nadj.21 = methods::as(nadj.21, "dgCMatrix")

      # degree-normalized transition matrix with jumping probability
      nadjM = rbind(cbind(nadj.1, nadj.12), cbind(nadj.21, nadj.2))
      # Ensuring columns sum to 1
      nadjM = sum2one(nadjM)
    }else if(norm == "degree"){
      # PPIs
      nadj.1 = matrix(0, nrow = nrow(adj.1), ncol = nrow(adj.1), dimnames = dimnames(adj.1))
      c.sum = Matrix::colSums(adj.1)
      jump.tmp = ifelse(bp.1, 1 - jump, 1)
      for(j in which(c.sum != 0)){
        id.tmp = adj.1@i[(adj.1@p[j] + 1):adj.1@p[j+1]] + 1 # row ids of non-zero elements in column j
        nadj.1[id.tmp, j] = jump.tmp[j] * adj.1[id.tmp, j] / c.sum[j]
      }
      nadj.1 = methods::as(nadj.1, "dgCMatrix")

      # MMIs
      nadj.2 = matrix(0, nrow = nrow(adj.2), ncol = nrow(adj.2), dimnames = dimnames(adj.2))
      c.sum = Matrix::colSums(adj.2)
      jump.tmp = ifelse(bp.2, 1 - jump, 1)
      for(j in which(c.sum != 0)){
        id.tmp = adj.2@i[(adj.2@p[j] + 1):adj.2@p[j+1]] + 1 # row ids of non-zero elements in column j
        nadj.2[id.tmp, j] = jump.tmp[j] * adj.2[id.tmp, j] / c.sum[j]
      }
      nadj.2 = methods::as(nadj.2, "dgCMatrix")

      # pm: proteins in rows, metabolites in columns
      nadj.12 = matrix(0, nrow = nrow(adj.12), ncol = ncol(adj.12), dimnames = dimnames(adj.12))
      c.sum = Matrix::colSums(adj.12)
      for(j in which(c.sum != 0)){
        id.tmp = adj.12@i[(adj.12@p[j] + 1):adj.12@p[j+1]] + 1 # row ids of non-zero elements in column j
        nadj.12[id.tmp, j] <- jump * adj.12[id.tmp, j] / c.sum[j]
      }
      nadj.12 = methods::as(nadj.12, "dgCMatrix")

      # mp: metabolites in rows, proteins in columns
      nadj.21 = matrix(0, nrow = nrow(adj.21), ncol = ncol(adj.21), dimnames = dimnames(adj.21))
      c.sum = Matrix::colSums(adj.21)
      for(j in which(c.sum != 0)){
        id.tmp = adj.21@i[(adj.21@p[j] + 1):adj.21@p[j+1]] + 1 # row ids of non-zero elements in column j
        nadj.21[id.tmp, j] <- jump * adj.21[id.tmp, j] / c.sum[j]
      }
      nadj.21 = methods::as(nadj.21, "dgCMatrix")

      # degree-normalized transition matrix with jumping probability
      nadjM = rbind(cbind(nadj.1, nadj.12), cbind(nadj.21, nadj.2))
      # Ensuring columns sum to 1
      nadjM = sum2one(nadjM)
    }else{
      nadjM <- adjM
    }
    return(nadjM)
  }
}

#' @title Scale a transition matrix to be approximately bistochastic.
#'
#' @description This function uses Iterative Proportional Fitting to modify the input left-stochastic matrix such that the row and column sums all equal 1. The diagonal element in each column is then redistributed evenly to the other non-zero elements of that column.
#'
#' @details
#' The purpose of this function is to mitigate degree bias that is present in PPI networks by directly manipulating the transition matrix to be used in RWR.
#'
#' The stationary distribution of a Markov chain can be viewed as a measure of centrality. These are the values at which each node has total in-flow equal to total out-flow. Any initial distribution will converge to this distribution in the long run.
#'
#' For transition matrices obtained by column-normalizing an adjacency matrix, the stationary probabilities are proportional to node degree. If the goal is to minimize the influence that degree has on propagation scores, while still recognizing that all else being equal a higher-degree node should have more weight than a lower-degree node,
#' we should aim to squeeze stationary probabilities towards some global mean, since the stationary probability is a good proxy for the amount of influence of degree. The metric that best captures this goal is the entropy of the stationary distribution.
#'
#' A bistochastic (doubly stochastic) matrix has a uniform stationary distribution, which represents a stationary distribution with maximum entropy.
#' There are matrix scaling methods that modify a matrix to conform to the row and column sums of some target matrix, while still being as similar as possible to the original matrix. Here we use Iterative Proportional Fitting (IPF) to arrive at an approximately bistochastic matrix from an input transition matrix.
#'
#' To ensure convergence of IPF, it is necessary to add self-loops, i.e., non-zero diagonal elements. Interestingly, these diagonal elements mostly vanish to zero except for very low-degree nodes.
#' After obtaining an approximately bistochastic matrix, it is desirable to remove the self-loops, since this was only added for convergence considerations and the original matrix had zeros along the diagonal.
#' To do this, the diagonal elements are evenly redistributed to the non-zero elements of the column, thus upsetting the bistochastic approximation but preserving column sums. However, the resulting left-stochastic matrix is still much closer to being bistochastic than originally, and this is reflected by an increase in the entropy of the associated stationary distribution.
#'
#' @param trans_mat A transition matrix
#'
#' @return A left stochastic transition matrix
#'
#' @examples
#' # Calculate the entropy of a vector
#' entropy = function(x){
#'   tmp = ifelse(round(x, 20) == 0, 0, -x * log(x))
#'   sum(tmp)
#' }
#' # Calculate the stationary distribution of the transition matrix of an irreducible markov chain
#' stationary.distr = function(x){
#'   e = Re(RSpectra::eigs(A = x, k = 1, which = "LM")$vectors[,1])
#'   e / sum(e)
#' }
#' d = igraph::degree(glut4_graph) # node degrees of the graph
#' # normalizing adjacency matrix to get transition matrix
#' adj_norm = AMEND::transition_matrix(glut4_graph, glut4_adjM, norm = "degree")
#' p1 = stationary.distr(adj_norm)
#' e1 = entropy(p1)
#' # Perform bistochastic scaling
#' # This aims to maximize entropy of stationary distribution of transition matrix
#' adj_norm = AMEND::bistochastic_scaling(trans_mat = adj_norm)
#' p2 = stationary.distr(adj_norm)
#' e2 = entropy(p2)
#'
#' # Compare before and after
#' e1 < e2
#' plot(d, p1, xlab = "Degree", ylab = "Stationary Distribution",
#'      main = paste0("Before (entropy=", round(e1, 3),")"), ylim = c(0, 0.01))
#' abline(a = 0, b = 1 / sum(d), xpd = FALSE)
#' plot(d, p2, xlab = "Degree", ylab = "Stationary Distribution",
#'      main = paste0("After (entropy=", round(e2, 3),")"), ylim = c(0, 0.01))
#' abline(a = 0, b = 1 / sum(d), xpd = FALSE)
#'
#' @export
bistochastic_scaling = function(trans_mat){
  get_diagonal = function(X){
    if(!"dgCMatrix" %in% class(X)) X = Matrix::Matrix(X, sparse = TRUE)
    d = numeric(ncol(X))
    for(i in 1:ncol(X)){
      id.tmp = (X@p[i] + 1):X@p[i+1] # ids of X@x that are non-zero and in col i
      row.ids = X@i[id.tmp] + 1 # row ids of non-zero elements in col i
      if(i %in% row.ids){ # if diagonal is non-zero
        d[i] = X@x[id.tmp[row.ids == i]]
      }else next
    }
    d
  }
  # Iterative Proportional Fitting
  ipf = function(X, e = 1e-6){
    # For starting transition matrix X, obtain B = PXQ, where B is bistochastic, P,Q are diagonal matrices, and B and X are as similar as possible (minimize relative entropy between B & X)
    # Adding self-loops to aid in convergence
    d = get_diagonal(X)
    d[d == 0] = e
    X = X + Matrix::Diagonal(n = nrow(X), x = d)

    stop_delta = 1e-6
    step = 1
    stop_step = 200
    q1 = rep(1, nrow(X)) # initialize the diagonal elements of Q
    p1 = rep(1, nrow(X)) # initialize the diagonal elements of P
    repeat{ # To generalize to any row/col sums, change 1 to variables
      # set p, given q (p are the diagonal elements of P)
      p2 = 1 / Matrix::rowSums(X %*% Matrix::Diagonal(x = q1))
      # set q, given p found above (q are the diagonal elements of Q)
      q2 = 1 / Matrix::colSums(Matrix::Diagonal(x = p2) %*% X)

      delta.p = all(abs(p2 - p1) <= stop_delta)
      delta.q = all(abs(q2 - q1) <= stop_delta)
      step = step + 1
      if((delta.p && delta.q) || step > stop_step) break
      q1 = q2
      p1 = p2
    }
    P = Matrix::Diagonal(x = p2)
    Q = Matrix::Diagonal(x = q2)
    B = P %*% X %*% Q
    return(list(B = B, p = p2, q = q2))
  }
  x = 1e-6
  if(!"dgCMatrix" %in% class(trans_mat)) trans_mat = Matrix::Matrix(trans_mat, sparse = TRUE)
  B.tmp = ipf(trans_mat, x)
  B = B.tmp$B
  b = get_diagonal(B)
  tmp.res = numeric(Matrix::nnzero(B))
  for(i in 1:ncol(B)){
    if(b[i] == 0) next
    id.tmp = (B@p[i] + 1):B@p[i+1] # ids of B@x that are non-zero and in col i
    row.ids = B@i[id.tmp] + 1 # row ids of non-zero elements in col i
    diag.id = id.tmp[row.ids == i]
    off.id = id.tmp[row.ids != i]
    # Evenly distribute to neighbors of node i. Preserves column sums
    degr = B@p[i+1] - B@p[i] - 1 # number of non-zero elements in col i i.e., degree of node i. Minus 1 b/c of self-loops
    tmp = b[i] / degr
    tmp.res[diag.id] = 0
    tmp.res[off.id] = B@x[off.id] + tmp
  }
  B@x = tmp.res
  return(B)
}

#' @title Random walk with restart (RWR) procedure
#'
#' @description `RWR()` implements the RWR procedure. This simulates random walkers on a graph, with a certain probability of returning to the seeds nodes.
#'
#' @param nadjM Normalized adjacency matrix
#' @param setSeeds Vector of seed values
#' @param restart Restart parameter
#' @param heterogeneous Logical. If TRUE, network is considered heterogeneous (two distinct node types, e.g., proteins and metabolites). If TRUE, node_type must be included as an argument or graph vertex attribute
#' @param node_type Named vector of node types. Only considered if heterogeneous=TRUE. Can only accommodate two unique node types.
#' @param net1.weight Relative weight to give to nodes of network 1 (the first type to appear in node_type).
#'
#' @return vector of propagation scores from RWR
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
#'
#' # Normalizing the adjacency matrix
#' adj_norm = adjm %*% diag(1 / degree(g))
#'
#' # Creating a named vector of seed values
#' seeds = runif(8)
#' names(seeds) = 1:8
#'
#' rwr = AMEND:::RWR(nadjM = adj_norm, setSeeds = seeds, restart = 0.8)
#'
RWR <- function (nadjM, setSeeds = NULL, restart = 0.75, heterogeneous = FALSE, node_type = NULL, net1.weight = 0.5){
  # NB: Assuming node_type in same order as setSeeds, while order of setSeeds doesn't need to equal order of rows/columns in nadjM
  # The ordering will conform to nadjM
  # Also, assuming that node types in nadjM are in blocks (e.g., all proteins, then all metabolites)

  sum2one <- function(X) {
    if(!"dgCMatrix" %in% class(X)) X = Matrix::Matrix(X, sparse = TRUE)

    inv_col_sum = Matrix::Diagonal(x = abs(Matrix::colSums(X))^(-1))
    res = X %*% inv_col_sum
    res[is.na(res)] = 0
    res
  }
  if (is.null(restart) || is.na(restart) || restart < 0 || restart > 100) {
    r <- 0.75
  }else if (restart > 1 && restart < 100) {
    r <- restart/100
  }else {
    r <- restart
  }
  stop_delta <- 1e-06
  stop_step <- 50
  if (is.null(setSeeds)) {
    P0matrix <- Matrix::Matrix(diag(nrow(nadjM)), sparse = T)
    rownames(P0matrix) <- rownames(nadjM)
    colnames(P0matrix) <- rownames(nadjM)
  }else {
    if(heterogeneous && is.null(node_type)) stop("node_type is not given. node_type must be given for heterogeneous graphs.")
    if (is.matrix(setSeeds) | is.data.frame(setSeeds)) {
      data <- as.matrix(setSeeds)
    }else if (is.vector(setSeeds)) {
      data <- as.matrix(setSeeds, ncol = 1)
    }
    if (is.null(rownames(data))) {
      stop("The function must require the row names of the input setSeeds.\n")
    }else if (any(is.na(rownames(data)))) {
      warning("setSeeds with NA as row names will be removed")
      data <- data[!is.na(rownames(data)), ]
      node_type = node_type[!is.na(rownames(data))]
    }
    ind <- match(rownames(data), rownames(nadjM)) # This ensures that all(rownames(data) == rownames(nadjM)[ind]) is TRUE
    nodes_mapped <- rownames(nadjM)[ind[!is.na(ind)]]
    node_type = node_type[!is.na(ind)]
    nt = unique(node_type)
    if (length(nodes_mapped) != nrow(nadjM)) {
      warning("The row names of input setSeeds do not contain all those in the input graph.\n")
    }

    P0matrix <- matrix(0, nrow = nrow(nadjM), ncol = ncol(data))
    P0matrix[ind[!is.na(ind)], ] <- as.matrix(data[!is.na(ind),]) # ensuring that order of P0matrix is same as rownames(nadjM)
    if(heterogeneous){
      node_type[ind[!is.na(ind)]] <- node_type # ensuring that order of node_type is same as P0matrix
      # Setting seed vector
      P0matrix <- rbind(net1.weight * sum2one(P0matrix[node_type == nt[1],]),
                        (1 - net1.weight) * sum2one(P0matrix[node_type == nt[2],]))
    }else{
      P0matrix <- sum2one(P0matrix)
    }
    # P0matrix <- Matrix::Matrix(P0matrix, sparse = T)
  }
  if (restart == 1) {
    PTmatrix <- P0matrix
  }else {
    PTmatrix <- Matrix::Matrix(0, nrow = nrow(P0matrix), ncol = ncol(P0matrix), sparse = T)
    P0 <- P0matrix[, 1]
    step <- 0
    delta <- 1
    PT <- P0
    while (delta > stop_delta && step <= stop_step) {
      PX <- (1 - r) * nadjM %*% PT + r * P0
      delta <- sum(abs(PX - PT))
      PT <- PX
      step <- step + 1
    }
    PT[PT < 1e-10] <- 0
    PTmatrix[, 1] <- PT
  }
  PTmatrix <- sum2one(PTmatrix)
  PTmatrix[PTmatrix < 1e-10] <- 0
  # PTmatrix <- Matrix::Matrix(PTmatrix, sparse = T)
  rownames(PTmatrix) <- rownames(P0matrix)
  colnames(PTmatrix) <- colnames(P0matrix)
  invisible(PTmatrix)
}

#' @title Heuristic solution of Maximum-weight Connected Subgraph problem
#'
#' @description Given a graph and a named vector of node scores, `heinz()` heuristically finds a solution to the maximum-weight connected subgraph problem.
#'
#' @details Details can be found in the _dnet_ package documentation \url{https://cran.r-project.org/web/packages/dnet/dnet.pdf}
#'
#' @note Adapted from the _dnet_ package. Based on method from the _BioNet_ package
#' @param g Input graph
#' @param scores Named vector of scores for nodes in g
#'
#' @return a subnetwork as an igraph object
#'
#' @examples
#' library(igraph)
#'
#' graph = sample_pa(n = 100, power = 1.2)
#' V(graph)$name = 1:100
#' g_scores = rnorm(100)
#' names(g_scores) = 1:100
#'
#' new_g = AMEND:::heinz(g = graph, scores = g_scores)
#'
heinz <- function (g, scores){
  ig <- g
  if(class(ig) != "igraph"){
    stop("The function must apply to an 'igraph' object.\n")
  }
  if(is.null(V(ig)$name)){
    V(ig)$name <- as.character(V(ig))
  }
  if(is.null(names(scores))){
    stop("The function must require the names of the input scores.\n")
  }else if(any(is.na(names(scores)))){
    warning("Those scores with NA as names will be removed")
    scores <- scores[!is.na(names(scores))]
  }
  V(ig)$score <- scores[V(ig)$name]
  pos.nodes = which(V(ig)$score > 0)
  if(length(pos.nodes) == 0){
    warning("No positive nodes")
    subgraph <- igraph::graph.empty(n = 0, directed = F)
  }else if(length(pos.nodes) == 1){
    subgraph = igraph::induced_subgraph(ig, pos.nodes)
    V(subgraph)$type <- "desired"
  }else{
    pos.subgraph = igraph::induced_subgraph(ig, pos.nodes)
    conn.comp.graph <- igraph::decompose(pos.subgraph)
    score.comp <- unlist(lapply(lapply(conn.comp.graph,
                                       function(x) as.numeric(V(x)$score)), sum))
    ind_order <- order(score.comp, decreasing = T)
    conn.comp.graph <- conn.comp.graph[ind_order]
    score.comp <- score.comp[ind_order]
    for(i in 1:length(conn.comp.graph)){
      conn.comp.graph[[i]]$score <- score.comp[i]
    }
    v.id <- seq(1, vcount(ig))
    names(v.id) <- V(ig)$name
    edgelist <- igraph::get.edgelist(ig, names = F)
    edgelist1 <- edgelist[, 1]
    edgelist2 <- edgelist[, 2]
    #===#
    # This for loop is changing the edgelist to treat the nodes in each meta-node as one node
    ig.size <- vcount(ig)
    for(i in 1:length(conn.comp.graph)){ # for each meta-node i
      # change the source, target ids of all edges containing nodes in meta-node i to be the same
      new.id <- ig.size + i
      v.id.tmp = v.id[V(conn.comp.graph[[i]])$name]
      edgelist1[edgelist1 %in% v.id.tmp] = new.id
      edgelist2[edgelist2 %in% v.id.tmp] = new.id
    }
    #===#
    new.ids <- seq(vcount(ig) + 1, vcount(ig) + length(conn.comp.graph))
    new.names <- paste("cluster", seq(1:length(conn.comp.graph)), sep = "")
    names(new.ids) <- new.names
    v.id <- c(v.id, new.ids)
    v.name <- names(v.id)
    names(v.name) <- v.id
    new.edgelist <- cbind(v.name[as.character(edgelist1)],
                          v.name[as.character(edgelist2)])
    mig <- igraph::graph_from_edgelist(new.edgelist, directed = F)
    mig <- igraph::simplify(mig, remove.loops = T, remove.multiple = T)
    node.score <- scores[V(mig)$name]
    names(node.score) <- V(mig)$name
    node.score.cluster <- sapply(conn.comp.graph, igraph::get.graph.attribute, "score")
    names(node.score.cluster) <- new.names
    ind_cluster <- grep("cluster", names(node.score))
    node.score[ind_cluster] <- node.score.cluster[names(node.score[ind_cluster])]
    V(mig)$score <- node.score
    score.degree <- 1/(igraph::degree(mig) + 1)
    tmp_score <- V(mig)$score
    tmp_score[tmp_score > 0] <- 0
    V(mig)$score.degree <- score.degree * tmp_score
    E(mig)$weight <- rep(0, length(E(mig)))
    tmp_n1 <- igraph::get.edgelist(mig, names = F)[, 1]
    tmp_n2 <- igraph::get.edgelist(mig, names = F)[, 2]
    E(mig)$weight <- -(V(mig)[tmp_n1]$score.degree + V(mig)[tmp_n2]$score.degree)
    if(!igraph::is_connected(mig)){
      decomp.graphs <- igraph::decompose(mig)
      sum.pos <- lapply(decomp.graphs, function(x) {
        sum(node.score[names(which(node.score[V(x)$name] > 0))])
      })
      mig <- decomp.graphs[[which.max(sum.pos)]]
      rm(decomp.graphs)
    }
    mst <- igraph::minimum.spanning.tree(mig, weights = E(mig)$weight)
    mst.cluster.id <- grep("cluster", V(mst)$name)
    names(mst.cluster.id) <- V(mst)[mst.cluster.id]$name
    tmp <- unlist(strsplit(names(mst.cluster.id), "cluster"))
    ttmp <- as.numeric(matrix(tmp, nrow = 2)[2, ])
    mst.cluster.id <- mst.cluster.id[order(ttmp)]
    if(length(mst.cluster.id) == 1){
      neg.node.ids.2 = c()
    }else{
      # Iteratively exclude all negative nodes of degree 1
      loop = TRUE
      while(loop){
        neg_deg1 = V(mst)[igraph::degree(mst) == 1 & V(mst)$score < 0]
        if(length(neg_deg1) == 0){
          loop = FALSE
        }else{
          mst = igraph::induced_subgraph(mst, -neg_deg1)
        }
      }
      sub.mig = igraph::induced_subgraph(mig, which(V(mig)$name %in% V(mst)$name))
      if(!igraph::is_simple(sub.mig)) sub.mig = igraph::simplify(sub.mig)
      neg.node.ids <- which(V(sub.mig)$score < 0)
      for(i in neg.node.ids){
        tmp_nei <- igraph::neighbors(sub.mig, v = i)
        tmp_nei_meta <- grep("cluster", V(sub.mig)[tmp_nei]$name)
        V(sub.mig)[i]$clusters <- list(tmp_nei[tmp_nei_meta])
      }
      score.neg.nodes <- c()
      for(i in neg.node.ids){
        if(!is.na(V(sub.mig)[i]$clusters[1])){
          borders <- c(i, V(sub.mig)[i]$clusters)
          borders <- unlist(borders)
          score.neg.nodes <- c(score.neg.nodes, sum(V(sub.mig)[borders]$score))
        }else{
          score.neg.nodes <- c(score.neg.nodes, V(sub.mig)[i]$score)
        }
      }
      neg.node.ids.2 <- neg.node.ids[score.neg.nodes > 0] # > or >= ?
    }
    if(length(neg.node.ids.2) == 0){
      tmp <- unlist(strsplit(names(node.score.cluster)[which.max(node.score.cluster)], "cluster"))
      ttmp <- as.numeric(matrix(tmp, nrow = 2)[2, ])
      tmp_nodes <- unlist(lapply(conn.comp.graph, igraph::get.vertex.attribute, "name")[ttmp])
      subgraph = igraph::induced_subgraph(ig, which(V(ig)$name %in% tmp_nodes))
      if(!igraph::is_connected(subgraph)){
        clust <- igraph::components(subgraph)
        cid <- which.max(clust$csize)
        subgraph <- igraph::induced_subgraph(subgraph, which(clust$membership == cid))
      }
      V(subgraph)$type <- "desired"
    }else{
      subg = igraph::induced_subgraph(sub.mig, neg.node.ids.2)
      if(!igraph::is_connected(subg)){
        clust <- igraph::components(subg)
        cid <- which.max(clust$csize)
        subg <- igraph::induced_subgraph(subg, which(clust$membership == cid))
      }
      if(!igraph::is_simple(subg)) subg = igraph::simplify(subg)
      mst.subg <- igraph::minimum.spanning.tree(subg, E(subg)$weight)
      getPathScore <- function(path, graph1, graph2) {
        s1 <- V(graph1)[path]$score
        tmp <- unique(unlist(V(graph1)[path]$clusters))
        s2 <- V(graph2)[tmp]$score
        sum(c(s1, s2))
      }
      max.score <- 0
      best.path <- c()
      for(i in 1:vcount(mst.subg)){
        path <- igraph::all_shortest_paths(mst.subg, from = V(mst.subg)[i])
        path.score <- unlist(lapply(path$res, getPathScore, graph1 = mst.subg, graph2 = sub.mig))
        best.pos <- which.max(path.score)
        if(path.score[[best.pos]] > max.score){
          best.path <- path$res[[best.pos]]
          max.score <- path.score[[best.pos]]
        }
      }
      if(length(best.path) != 1){
        cluster.list <- V(mst.subg)[best.path]$clusters
        names.list <- as.character(1:length(cluster.list))
        names(cluster.list) <- names.list
        names(best.path) <- names.list
        for(i in names.list){
          res <- lapply(cluster.list, intersect, cluster.list[[i]])
          if(length(intersect(unlist(cluster.list[as.character(which(as.numeric(names.list) < as.numeric(i)))]),
                              unlist(cluster.list[as.character(which(as.numeric(names.list) > as.numeric(i)))]))) > 0){
            if(length(setdiff(res[[i]], unique(unlist(res[names(res) != i])))) == 0){
              cluster.list <- cluster.list[names(cluster.list) != i]
              names.list <- names.list[names.list != i]
            }
          }
        }
        best.path <- best.path[names.list]
      }
      pos.cluster <- V(sub.mig)[unique(unlist(V(mst.subg)[best.path]$clusters))]$name
      tmp <- unlist(strsplit(pos.cluster, "cluster"))
      ttmp <- as.numeric(matrix(tmp, nrow = 2)[2, ])
      tmp_meta_nodes <- unlist(lapply(conn.comp.graph, igraph::get.vertex.attribute, "name")[ttmp])
      tmp_border_nodes <- V(mst.subg)[best.path]$name
      tmp_nodes <- c(tmp_border_nodes, tmp_meta_nodes)
      subgraph <- igraph::induced_subgraph(ig, vids = tmp_nodes)
      if(!igraph::is_connected(subgraph)){
        clust <- igraph::components(subgraph)
        cid <- which.max(clust$csize)
        subgraph <- igraph::induced_subgraph(subgraph, which(clust$membership == cid))
      }
      type <- rep("desired", vcount(subgraph))
      names(type) <- V(subgraph)$name
      type[tmp_border_nodes[tmp_border_nodes %in% names(type)]] <- "linker"
      V(subgraph)$type <- type
    }
  }
  return(subgraph)
}

#' @title Exponential decay schedule for filtering rate
#'
#' @description Exponential decay schedule for filtering rate. Used to shift the raw RWR scores
#'
#' @param eta0 Starting filtering rate
#' @param d Decay parameter
#'
#' @return vector of filtering rates for each iteration of AMEND
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
#' @return decay value, numeric
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
#' @return Vector of core-clustering coefficients for each node
#'
#' @examples
#' core_clust_coefs = AMEND:::core_cc(glut4_graph)
#' head(core_clust_coefs)
#'
core_cc = function(g, weight = TRUE){
  # if(!igraph::is_connected(g)) stop("Input graph must be connected")
  if(!"weight" %in% igraph::edge_attr_names(g) && weight){
    # warning("No \'weight\' edge attribute. Treating as unweighted.")
    weight = FALSE
  }
  n = vcount(g)
  w = numeric(n)
  for(i in seq_len(n)){
    nborhood = igraph::make_ego_graph(g, order = 1, nodes = i)[[1]]
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
#' @return weighted edge density of graph
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
    N = igraph::vcount(graph)
    return(sum(igraph::E(graph)$weight) / (N * (N - 1) / 2))
  }else return(igraph::edge_density(graph))
}

#' @title Get a larger subnetwork from previous iterations of AMEND
#'
#' @description
#'
#' The AMEND algorithm implemented in 'run_AMEND()' may overshoot the approximate final module size, or the user may be interested in the parent subnetworks of the final module. This function allows the user to examine intermediate subnetworks that were generated during 'run_AMEND()'.
#'
#' @param amend_object Result from `run_AMEND()`
#' @param k Iteration associated with the subnetwork you want to retrieve
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
#' subnet1 = run_AMEND(graph = glut4_graph, data.type = "ECI", DOI = "negative")
#'
#' # Inspect the statistics of the intermediate subnetworks
#' subnet1$stats
#'
#' # Retrieve the subnetwork of iteration 7
#' subnet7 = get_subnetwork(subnet1, 7)
#' subnet7
#'
#' @export
get_subnetwork = function(amend_object, k){
  ig = amend_object$input_params$ig
  nodes = amend_object$subnetworks[[k]]
  igraph::induced_subgraph(ig, which(V(ig)$name %in% nodes))
}
