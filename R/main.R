# Links for R package development
# https://r-pkgs.org/
# https://kbroman.org/pkg_primer/

# Link for Github guide
# https://kbroman.org/github_tutorial/


#' @importFrom igraph vcount V V<- E E<- vertex_attr vertex_attr<-
#' @importFrom stats sd quantile

#' @title Pipeline for AMEND workflow
#'
#' @description `run_AMEND()` performs all steps of the AMEND workflow, including setting the starting filtering rate through PSO and identifying an 'active' module.
#'
#' @details AMEND is an active module identification method. An active module is a subset of connected nodes in an interaction network that have relatively large experimental values.
#' AMEND takes as input a protein-protein interaction (PPI) network and node-wise summaries of experimental data (e.g., log fold change from a RNA-seq differential expression analysis). It returns a single connected subnetwork (i.e., module).
#'
#' The final module is found in an iterative manner by applying random walk with restart (RWR) on the network to get node weights, which are then shifted by a certain quantile of the RWR weights (termed the filtering rate) so that there are both negative and positive values.
#' These shifted RWR weights are used to find the maximum-weight connected subgraph. This new subgraph is input into RWR again and the process repeats until there is no change between iterations. The subnetwork with largest score is returned as the final module.
#'
#' The filtering rate changes at each iteration following an exponential decay schedule, with the starting filtering rate determined by particle swarm optimization (PSO).
#'
#' The network score is the product of the average standardized experimental values and the average core-clustering coefficient of the nodes in the network. The experimental values (e.g., log fold changes) are standardized w.r.t. the entire dataset.
#'
#' See \url{/link/to/paper}.
#'
#' @param graph an igraph object with experimental data as vertex attribute. Either graph or adj_matrix and node_scores must be specified
#' @param n (integer): approximate size of the final module. Used for setting the filtering rate schedule
#' @param eta (numeric): starting filtering rate. If NULL (default), determines optimal value through PSO
#' @param data.type (character): one of ECI, logFC, or p_val. If graph is specified, must have matching vertex attribute
#' @param eci.direction (character): direction of interest for ECI values. negative or positive
#' @param logFC.direction (character): direction of interest for log fold change values. negative, positive, or both. When "both" is specified, the absolute value of logFC is taken
#' @param adj_matrix (matrix): adjacency matrix of network. Either graph or adj_matrix and node_scores must be specified
#' @param node_scores vector of node scores
#' @param normalize (character): normalization scheme of adjacency matrix for random walk with restart
#' @param seed.weight (numeric): Relative weight to give to nodes not in the direction of interest in random walk with restart, between 0 and 1
#' @param max_it (integer): maximum number of iterations for PSO
#' @param n_particles (integer): number of particles for PSO
#' @param random_seed (integer): random seed. If NULL (default), no seed is set by the user
#' @param verbose (logical): Whether to output current iteration number to show progress
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
#' # Use run_AMEND() with an igraph object
#' subnet1 = run_AMEND(graph = glut4_graph, data.type = "ECI", eci.direction = "negative",
#'                     max_it = 1, n_particles = 1, verbose = TRUE)
#'
#' # Use run_AMEND() with an adjacency matrix and a vector of node scores
#' subnet2 = run_AMEND(adj_matrix = glut4_adjM, node_scores = eci_scores, data.type = "ECI",
#'                     eci.direction = "negative", max_it = 1, n_particles = 1, verbose = TRUE)
#' }
#'
#' @export
run_AMEND = function(graph, n = 25, eta = NULL, data.type = c("ECI", "logFC", "p_val"),
                          eci.direction = c("positive", "negative"), logFC.direction = c("positive", "negative", "both"),
                          adj_matrix = NULL, node_scores = NULL, normalize = c("core", "degree"), seed.weight = 0.5,
                          max_it = 3, n_particles = 3, random_seed = NULL, verbose = F){
  data.type <- match.arg(data.type)
  eci.direction <- match.arg(eci.direction)
  logFC.direction <- match.arg(logFC.direction)
  normalize <- match.arg(normalize)

  if(is.null(eta)){
    if(!is.null(random_seed) & is.numeric(random_seed)) set.seed(random_seed)

    eta.search = pso::psoptim(par = NA, fn = amend,  eta.search = TRUE, graph = graph, n = n, data.type = data.type,
                              eci.direction = eci.direction, logFC.direction = logFC.direction,
                              adj_matrix = adj_matrix, node_scores = node_scores, normalize = normalize, seed.weight = seed.weight, verbose = verbose,
                              lower = 0.2, upper = 0.8,
                              control = list(trace = 1, fnscale = -1, maxit = max_it, s = n_particles, w = 0.729, c.p = 1.49445, c.g = 1.49445, v.max = 0.5, type = "SPSO2011"))
    e = eta.search$par
  }else e = eta

  # Running AMEND with final starting filtering rate (eta)
  subnet = amend(eta = e, graph = graph, n = n, data.type = data.type,
                 eci.direction = eci.direction, logFC.direction = logFC.direction, adj_matrix = adj_matrix,
                 node_scores = node_scores, normalize = normalize, seed.weight = seed.weight, eta.search = FALSE, verbose = verbose)
  return(subnet)
}
# EDIT 11/1: added examples, added details
# EDIT 11/4: changed random_seed argument and implementation


#' @title Identify active modules from an interaction network
#'
#' @description Identifies an active module, i.e., subnetwork with large experimental values, by iteratively applying random walk with restart and a heurstic solution to the maximum-weight connected subgraph problem.
#'
#' @details See `run_AMEND()` and \url{link/to/paper}.
#'
#' @param eta (numeric): starting filtering rate
#' @param graph an igraph object with experimental data as vertex attribute. Either graph or adj_matrix and node_scores must be specified
#' @param n (integer): approximate size of the final module
#' @param data.type (character): one of ECI, logFC, or p_val. If graph is specified, must have matching vertex attribute
# @param data_prep_func (function): function for preprocessing data. First argument should be an igraph object, it must output an igraph object with vertex attributes "seeds" and "Z", and should allow for optional arguments
#' @param eci.direction (character): direction of interest for ECI values. negative or positive
#' @param logFC.direction (character): direction of interest for log fold change values. negative, positive, or both. When "both" is specified, the absolute value of logFC is taken
#' @param adj_matrix (matrix): adjacency matrix of network. Either graph or adj_matrix and node_scores must be specified
#' @param node_scores vector of node scores
#' @param normalize (character): normalization scheme of adjacency matrix for random walk with restart
#' @param eta.search (logical): Whether the starting filtering rate (eta) is being searched for
#' @param seed.weight (numeric): Relative weight to give to nodes not in the direction of interest in random walk with restart, between 0 and 1
#' @param verbose (logical): Whether to output current iteration number to show progress
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
#' # Use amend() with an igraph object with a vertex attribute matching data.type arg
#' subnet1 = amend(graph = glut4_graph, eta = 0.6, data.type = "ECI",
#'                 eci.direction = "negative", verbose = TRUE)
#'
#' # Use amend() with an adjacency matrix and a vector of node scores
#' subnet2 = amend(adj_matrix = glut4_adjM, node_scores = eci_scores, eta = 0.6,
#'                 data.type = "ECI", eci.direction = "negative", verbose = TRUE)
#' }
#'
#' @export
amend <- function(eta, graph, n = 25, data.type = c("ECI", "logFC", "p_val"),
                  eci.direction = c("positive", "negative"), logFC.direction = c("positive", "negative", "both"),
                  adj_matrix = NULL, node_scores = NULL, normalize = c("core", "degree"), eta.search = FALSE,
                  seed.weight = 0.5, verbose = F){
  start_time <- Sys.time()

  data.type <- match.arg(data.type)
  eci.direction <- match.arg(eci.direction)
  logFC.direction <- match.arg(logFC.direction)
  normalize <- match.arg(normalize)

  if(all(class(adj_matrix) != "matrix") & !is.null(adj_matrix)){
    adj_matrix <- as.matrix(adj_matrix)
  }

  if(!is.null(adj_matrix) & !is.null(node_scores)){
    graph = igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = "weight")
    graph = igraph::set_vertex_attr(graph = graph, name = data.type, value = node_scores)
  }else if(!data.type %in% igraph::vertex_attr_names(graph)){
    stop(paste0("Graph has no vertex attribute \"", data.type, "\""))
  }

  # Determining if graph has edge attribute 'weight'
  if(!"weight" %in% igraph::edge_attr_names(graph)) E(graph)$weight = 1

  # Transforming input data
  graph = data_preprocessing(graph, data.type, eci.direction, logFC.direction, seed.weight)
  # data_prep_fun = match.fun(data_prep_fun)
  # graph = data_prep_fun(graph, data.type, eci.direction, logFC.direction, seed.weight, ...)

  repeat{
    if(eta.search) message(paste('Starting filtering rate:', round(eta, 3)))

    all_scores <- list()
    all_nets <- list()

    sn_score <- c(-Inf, 0)
    subg <- vector(mode = "list", length = 2)
    subg[[1]] <- graph
    k <- 1

    decay <- find.decay(vcount(graph), eta0 = eta, n = n)
    l <- exp.filtering.rate(eta0 = eta, d = decay)

    repeat{
      if(verbose) message(paste("Iteration:", k))
      if(k == 1){
        if(is.null(adj_matrix)){
          adj_og <- igraph::as_adjacency_matrix(subg[[1]], type = "both", attr = "weight", sparse = T)
          adj <- adj_og
        }else{
          adj_og <- adj_matrix
          adj <- adj_og
        }
      }else{
        adj <- adj[-rm.id, -rm.id]
      }

      # Normalize adjacency matrix
      n.adjM <- norm.adj(g = subg[[1]], adjM = adj, norm = normalize)

      Seeds <- data.frame(Experimental.values = V(subg[[1]])$seeds, row.names = V(subg[[1]])$name)

      # Choosing Restart parameter value through multi-level grid search
      # The search space changes depending on size of network
      if(k == 1){
        r = 0.99

        PTmatrix2 <- dRWR_alt(g = subg[[1]], nadjM = n.adjM, setSeeds = Seeds, restart = r)
        rwr_res <- as.vector(PTmatrix2)

        # Scores for Maximum Scoring Subgraph Algorithm
        rwr_score <- rwr_res - quantile(rwr_res, l[k])
        score <- rwr_score
        names(score) <- V(subg[[1]])$name

        subg[[2]] <- dNetFind_alt(subg[[1]], score)
        n.v <- vcount(subg[[2]])
        mZ <- sum(V(subg[[2]])$Z) / n.v
        mCC = sum(core_cc(subg[[2]])) / n.v
        score_metrics <- matrix(c(mCC, mZ, n.v), nrow = 1, ncol = 3)
        max.id <- 1
      }else if(vcount(subg[[1]]) >= 1000){
        srgs <- smart.restart.grid.search(0.7, 0.99, 10, 3, subg[[1]], n.adjM, Seeds, l[k])
        subg[[2]] <- srgs[[1]]
        score_metrics <- matrix(srgs[[2]], nrow = 1, ncol = 3)
        r <- srgs[[3]]
        max.id <- 1
      }else{
        srgs <- smart.restart.grid.search(0.2, 0.99, 18, 3, subg[[1]], n.adjM, Seeds, l[k])
        subg[[2]] <- srgs[[1]]
        score_metrics <- matrix(srgs[[2]], nrow = 1, ncol = 3)
        r <- srgs[[3]]
        max.id <- 1
      }
      # break out of repeat loop if there is no change in network
      if(vcount(subg[[2]]) == vcount(subg[[1]]) | vcount(subg[[2]]) <= 2) break

      mCC <- score_metrics[max.id, 1]
      mZ <- score_metrics[max.id, 2]
      n_nodes = score_metrics[max.id, 3]

      sn_score[2] <- mCC * mZ

      all_scores[[length(all_scores) + 1]] <- c(decay, r, sn_score[2], mZ, mCC, n_nodes, igraph::ecount(subg[[2]]), igraph::edge_density(subg[[2]]), l[k], mean(!V(subg[[1]])$name %in% V(subg[[2]])$name), seed.weight)
      all_nets[[length(all_nets) + 1]] <- subg[[2]]

      rm.id <- which(!V(subg[[1]])$name %in% V(subg[[2]])$name)

      k <- k + 1
      subg[[1]] <- subg[[2]]
      sn_score[1] <- sn_score[2]

    }
    message("*** Converged! ***")

    allscores <- do.call("rbind", all_scores)
    colnames(allscores) <- c("Decay", "Restart parameter", "Network score", "Avg Z", "Avg CCC", "Nodes", "Edges", "Density",
                             "Filtering rate", "Observed filtering rate", "Seed weight")
    allscores <- round(as.data.frame(allscores), 3)

    break
  }
  # Choosing subnetwork with maximum score
  net_ids = lapply(all_nets, function(x){
    row.names(adj_og)[row.names(adj_og) %in% V(x)$name]
    })

  max.id <- which.max(allscores[,3])
  best_sn <- all_nets[[max.id]]
  best_score <- allscores[max.id,3]

  end_time <- Sys.time()
  time <- end_time - start_time

  rm(list = ls()[!ls() %in% c("eta", "graph", "data.type", "eci.direction", "logFC.direction", "adj_matrix", "normalize",
                              "eta.search", "seed.weight", "best_sn", "best_score", "allscores", "time", "net_ids")])

  input_params = list(eta = eta, ig = graph, data.type = data.type, eci.dir = eci.direction, logFC.dir = logFC.direction, adjM = adj_matrix,
                      norm = normalize, eta.search = eta.search, seed.weight = seed.weight)

  if(!eta.search){
    return(list(module = best_sn, score = best_score, subnetworks = net_ids, stats = allscores, time = time, input_params = input_params))
  }else{
    rm(list = ls()[!ls() %in% c("best_score")])
    return(best_score)
  }
}
# EDIT 11/1: added examples, changed "weighted" arg in graph_from_adjacency_matrix(), modified to accomodate unweighted graphs, modified details
# EDIT 11/7: changed stopping criteria

#' @title Pre-processing of experimental data for input into AMEND and RWR
#'
#' @description Pre-processing of experimental data for input into AMEND and RWR. For use inside `amend()` and `run_AMEND()`.
#'
#' @details Transforms experimental data for use with Random walk with restart (RWR), which requires non-negative seed values with at least one non-zero element. For subnetwork scoring, the experimental scores need to be standardized w.r.t. all nodes in the input network.
#' For data.type "ECI", the seeding scheme involves taking the absolute value of the ECIs then weighting those nodes NOT in the direction of interest by 'seed.weight'
#' For data.type "logFC", the log fold changes are multiplied by -1 if the direction of interest is negative, 1 if the direction of interest is positive, and the absolute value is taken if the direction of interest is both. Then the values are exponentiated
#' For data.type "p_val", the p-values are transformed by -log10(p-value)
#'
#' @param graph (igraph object): input graph
#' @param data.type (character): one of ECI, logFC, or p_val. If graph is specified, must have matching vertex attribute
#' @param eci.direction (character): direction of interest for ECI values. negative or positive
#' @param logFC.direction (character): direction of interest for log fold change values. negative, positive, or both. When "both" is specified, the absolute value of logFC is taken
#' @param seed.weight (numeric): Relative weight to give to nodes not in the direction of interest in random walk with restart, between 0 and 1
#'
#' @return igraph object with 'seeds' and 'Z' vertex attributes
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
#'                                         eci.direction = "negative", seed.weight = 0.5)
#' new_graph2 = AMEND:::data_preprocessing(graph = glut4_graph, data.type = "ECI",
#'                                         eci.direction = "negative", seed.weight = 0.2)
#'
#' # Comparing
#' plot(density(V(new_graph1)$seeds))
#' plot(density(V(new_graph2)$seeds))
#'
#' plot(V(new_graph1)$ECI, V(new_graph1)$Z, type = "b",
#'      xlim = c(-1, 1), xlab = "ECI", ylab = "Standardized ECI")
#'
data_preprocessing = function(graph, data.type = c("ECI", "logFC", "p_val"), eci.direction = c("positive", "negative"),
                              logFC.direction = c("positive", "negative", "both"), seed.weight = 0.5){
  data.type <- match.arg(data.type)
  eci.direction <- match.arg(eci.direction)
  logFC.direction <- match.arg(logFC.direction)

  # Transforming input data
  if(data.type == "ECI"){
    if(seed.weight < 0){
      stop("Seed weight needs to be non-negative. Should be between 0 and 1.")
    }else if(seed.weight > 1){
      warning("Seed weight should be between 0 and 1.")
    }
    s <- ifelse(eci.direction == "positive", 1, -1)
    vertex_attr(graph, "seeds") = ifelse(s == rep(1, vcount(graph)), ifelse(V(graph)$ECI > 0, V(graph)$ECI, seed.weight * abs(V(graph)$ECI)), ifelse(V(graph)$ECI < 0, abs(V(graph)$ECI), seed.weight * V(graph)$ECI))
    # vertex_attr(graph, "seeds") <- s*V(graph)$ECI + 1
    V(graph)$Z <- (s*V(graph)$ECI - mean(s*V(graph)$ECI))/sd(s*V(graph)$ECI)
  }else if(data.type == "logFC"){
    if(logFC.direction == "both"){
      vertex_attr(graph, "seeds") <- exp(abs(V(graph)$logFC))
      V(graph)$Z <- (abs(V(graph)$logFC) - mean(abs(V(graph)$logFC)))/sd(abs(V(graph)$logFC))
    }else if(logFC.direction %in% c("positive", "negative")){
      s <- ifelse(logFC.direction == "positive", 1, -1)
      vertex_attr(graph, "seeds") <- exp(s*V(graph)$logFC)
      V(graph)$Z <- (s*V(graph)$logFC - mean(s*V(graph)$logFC))/sd(s*V(graph)$logFC)
    }
  }else if(data.type == "p_val"){
    vertex_attr(graph, "seeds") <- -log(vertex_attr(graph, "p_val") + 0.001, base = 10)
    vertex_attr(graph, "Z") <- (V(graph)$seeds - mean(V(graph)$seeds))/sd(V(graph)$seeds)
  }
  return(graph)
}
# EDIT 11/1: changed log1p to log(x + 0.001, base = 10), added examples


#' @title Grid search for the restart parameter in RWR
#'
#' @description Grid search for the restart parameter in RWR. For use inside `amend()` and `run_AMEND()`
#'
#' @details This is a multi-level grid search. The first level is simply a grid search. For the subsequent levels, the neighborhood about the optimal grid value from the previous level is searched. Stops when no improvement in network score or maximum number of levels is reached.
#'
#' @param g.min (numeric): minimum grid value
#' @param g.max (numeric): maximum grid value
#' @param n.1 (integer): number of grid points to search at first level
#' @param levels (integer): number of levels to search
#' @param ig input graph
#' @param n.adj.M (matrix): normalized adjacency matrix
#' @param seeds vector of seed values
#' @param filtering_rate (numeric): quantile for shifting the raw RWR scores
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
#'                0, 0, 0, 0, 1, 0, 0, 0), nrow = 8)
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
#' search = AMEND:::smart.restart.grid.search(g.min = 0.1, g.max = 0.9, n.1 = 10, levels = 3,
#'                                            ig = g, n.adj.M = adj_norm, seeds = seeds,
#'                                            filtering_rate = 0.3)
#'
smart.restart.grid.search <- function(g.min, g.max, n.1, levels = 3, ig, n.adj.M, seeds, filtering_rate){
  # Create grid
  get.r.grid<- function(best.g, g.min, g.max, p1, n){
    p <- p1/(2^(n-1))
    if(best.g == g.min){
      new.grid <- best.g + p
    }else if(best.g == g.max){
      new.grid <- best.g - p
    }else new.grid <- c(best.g - p, best.g + p)
    return(new.grid)
  }

  p1 <- (g.max - g.min)/(n.1 - 1)
  grid <- seq(g.max, g.min, by = -p1)

  best.r <- numeric()
  best.scores <- c(-Inf)
  n <- 1

  all_metrics = c()

  repeat{
    score_metrics <- matrix(nrow = length(grid), ncol = 3)
    nets <- vector(mode = "list", length = length(grid))

    for(i in seq_along(grid)){
      # RWR
      PTmatrix2 <- dRWR_alt(g = ig, nadjM = n.adj.M, setSeeds = seeds, restart = grid[i])
      rwr_res <- as.vector(PTmatrix2)

      # Scores for Maximum Scoring Subgraph Algorithm
      rwr_score <- rwr_res - quantile(rwr_res, filtering_rate)
      score <- rwr_score
      names(score) <- V(ig)$name

      nets[[i]] <- dNetFind_alt(ig, score)

      n.v <- vcount(nets[[i]])

      mZ <- sum(V(nets[[i]])$Z) / n.v

      mCC = sum(core_cc(nets[[i]])) / n.v

      score_metrics[i,] <- c(mCC, mZ, n.v)
    }

    score_metrics = score_metrics[which(score_metrics[,1] != 0),]

    if(isTRUE(nrow(score_metrics) >= 1)){ # if score matrix isn't empty
      max.id <- which.max(score_metrics[, 2] * score_metrics[, 1])
      g.score <- score_metrics[max.id, 2] * score_metrics[max.id, 1]

      all_metrics = c(all_metrics, score_metrics[,3])

    }else if(n == 1){ # else if score matrix is empty on first level, break out of repeat loop
      new.g = ig
      n.v = vcount(new.g)
      new.scores = c(sum(core_cc(new.g)) / n.v, sum(V(new.g)$Z) / n.v, n.v)
      best.r = NA
      break
    }else g.score = -Inf # else if score matrix is empty NOT on first level

    if(g.score > best.scores[ifelse(n == 1, 1, n-1)]){ # if best of current level is better than previous level
      new.g <- nets[[max.id]]
      new.scores <- score_metrics[max.id,]
      best.scores[n] <- g.score
      best.r[n] <- grid[max.id]
    }else{ # else keep score from last level
      best.scores[n] <- best.scores[n-1]
      best.r[n] <- best.r[n-1]
    }

    if(n == levels) break

    n <- n + 1

    grid <- get.r.grid(best.r[n-1], g.min, g.max, p1, n)
  }

  return(list(new.g, new.scores, best.r[n]))
}
# EDIT 1/11: added examples, added details


#' @title Nnormalize the adjacency matrix of a graph
#'
#' @description `norm.adj()` normalizes the adjacency matrix of a graph for input into random walk with restart (RWR).
#'
#' @details There are two normalization schemes available.
#' * core: uses the concept of node coreness. Details can be found in \url{https://academic.oup.com/nar/article/48/17/e98/5879427}
#' * degree: divides each column of the adjacency matrix by the column total
#'
#' @param g (igraph object): input graph
#' @param adjM (matrix): adjacency matrix
#' @param norm (character): normalization method
#'
#' @return normalized adjacency matrix
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
#'                0, 0, 0, 0, 1, 0, 0, 0), nrow = 8)
#' g = graph_from_adjacency_matrix(adjm, mode = 'undirected')
#'
#' adj_norm = AMEND:::norm.adj(g, adjm)
#'
norm.adj <- function(g, adjM, norm = c("core", "degree")){
  norm <- match.arg(norm)
  ig <- g

  if (norm == "core") {
    adjE <- igraph::as_adjacency_matrix(ig, type = "both", sparse = FALSE)
    core <- igraph::coreness(ig)
    nadjM <- matrix(0, nrow = nrow(adjE), ncol = nrow(adjE))
    for(j in 1:nrow(nadjM)){
      cond <- adjE[,j] != 0
      denom <- sum(core[cond])
      nadjM[cond, j] <- core[cond]/denom
    }
  }
  else if (norm == "degree") {
    wt <- Matrix::Diagonal(x = Matrix::colSums(adjM)^(-1))
    nadjM <- adjM %*% wt
  }
  else {
    nadjM <- adjM
  }
  return(nadjM)
}
# EDIT 11/1: removed edges arg in as_adjacency_matrix(), added examples, removed "weight" option in norm arg, modified degree normalization


#' @title Random walk with restart (RWR) procedure
#'
#' @description `dRWR_alt()` implements the RWR procedure. This simulates random walkers on a graph, with a certain probability of returning to the seeds nodes.
#'
#' @note Adapted from the _dnet_ package
#' @param g (igraph object): input graph
#' @param nadjM (matrix): normalized adjacency matrix
#' @param setSeeds (vector): vector of seed values
#' @param restart (numeric): restart parameter
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
#'                0, 0, 0, 0, 1, 0, 0, 0), nrow = 8)
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
#' rwr = AMEND:::dRWR_alt(g, adj_norm, seeds, 0.8)
#'
dRWR_alt <- function (g, nadjM, setSeeds = NULL, restart = 0.75)
{
  ig <- g

  sum2one <- function(PTmatrix) {
    col_sum <- abs(apply(PTmatrix, 2, sum))
    col_sum_matrix <- matrix(rep(col_sum, nrow(PTmatrix)),
                             ncol = ncol(PTmatrix), nrow = nrow(PTmatrix), byrow = T)
    res <- as.matrix(PTmatrix)/col_sum_matrix
    res[is.na(res)] <- 0
    return(res)
  }
  if (is.null(restart) || is.na(restart) || restart < 0 || restart > 100) {
    r <- 0.75
  }
  else if (restart > 1 && restart < 100) {
    r <- restart/100
  }
  else {
    r <- restart
  }
  stop_delta <- 1e-06
  stop_step <- 50
  if (is.null(setSeeds)) {
    P0matrix <- Matrix::Matrix(diag(vcount(ig)), sparse = T)
    rownames(P0matrix) <- V(ig)$name
    colnames(P0matrix) <- V(ig)$name
  }
  else {
    if (is.matrix(setSeeds) | is.data.frame(setSeeds)) {
      data <- as.matrix(setSeeds)
    }
    else if (is.vector(setSeeds)) {
      data <- as.matrix(setSeeds, ncol = 1)
    }
    if (is.null(rownames(data))) {
      stop("The function must require the row names of the input setSeeds.\n")
    }
    else if (any(is.na(rownames(data)))) {
      warning("setSeeds with NA as row names will be removed")
      data <- data[!is.na(rownames(data)), ]
    }
    cnames <- colnames(data)
    if (is.null(cnames)) {
      cnames <- seq(1, ncol(data))
    }
    ind <- match(rownames(data), V(ig)$name)
    nodes_mapped <- V(ig)$name[ind[!is.na(ind)]]
    if (length(nodes_mapped) != vcount(ig)) {
      warning("The row names of input setSeeds do not contain all those in the input graph.\n")
    }
    P0matrix <- matrix(0, nrow = nrow(nadjM), ncol = ncol(data))
    P0matrix[ind[!is.na(ind)], ] <- as.matrix(data[!is.na(ind),
                                                   ])
    P0matrix <- sum2one(P0matrix)
    P0matrix <- Matrix::Matrix(P0matrix, sparse = T)
  }
  if (restart == 1) {
    PTmatrix <- P0matrix
  }
  else {
    PTmatrix <- Matrix::Matrix(0, nrow = nrow(P0matrix),
                               ncol = ncol(P0matrix), sparse = T)
    for (j in 1:ncol(P0matrix)) {
      P0 <- P0matrix[, j]
      step <- 0
      delta <- 1
      PT <- P0
      while (delta > stop_delta && step <= stop_step) {
        PX <- (1 - r) * nadjM %*% PT + r * P0
        delta <- sum(abs(PX - PT))
        PT <- PX
        step <- step + 1
      }
      PT[PT < 1e-06] <- 0
      PTmatrix[, j] <- PT
    }
  }
  PTmatrix <- sum2one(PTmatrix)
  PTmatrix[PTmatrix < 1e-06] <- 0
  PTmatrix <- Matrix::Matrix(PTmatrix, sparse = T)
  rownames(PTmatrix) <- rownames(P0matrix)
  colnames(PTmatrix) <- colnames(P0matrix)
  invisible(PTmatrix)
}
# EDIT 1/11: added examples


#' @title Heuristic solution of Maximum-weight Connected Subgraph problem
#'
#' @description Given a graph and a named vector of node scores, `dNetFind_alt()` heuristically finds a solution to the maximum-weight connected subgraph problem.
#'
#' @details Details can be found in the _dnet_ package documentation \url{https://cran.r-project.org/web/packages/dnet/dnet.pdf}
#'
#' @note Adapted from the _dnet_ package. Based on method from the _BioNet_ package
#' @param g (igraph object): input graph
#' @param scores (vector): named vector of scores for nodes in g
#'
#' @return subnetwork as an igraph object
#'
#' @examples
#' library(igraph)
#'
#' graph = sample_pa(n = 100, power = 1.2)
#' V(graph)$name = 1:100
#' g_scores = rnorm(100)
#' names(g_scores) = 1:100
#'
#' new_g = AMEND:::dNetFind_alt(g = graph, scores = g_scores)
#'
dNetFind_alt <- function (g, scores)
{
  ig <- g
  if (class(ig) != "igraph") {
    stop("The function must apply to an 'igraph' object.\n")
  }
  if (is.null(V(ig)$name)) {
    V(ig)$name <- as.character(V(ig))
  }
  if (is.null(names(scores))) {
    stop("The function must require the names of the input scores.\n")
  }
  else if (any(is.na(names(scores)))) {
    warning("Those scores with NA as names will be removed")
    scores <- scores[!is.na(names(scores))]
  }
  V(ig)$score <- scores[V(ig)$name]
  pos.nodes <- V(ig)[V(ig)$score > 0]$name
  if (length(pos.nodes) == 0) {
    warning("No positive nodes")
    subgraph <- igraph::graph.empty(n = 0, directed = F)
  }
  else if (length(pos.nodes) == 1) {
    subgraph <- dnet::dNetInduce(ig, pos.nodes, knn = 0, remove.loops = T,
                           largest.comp = T)
    V(subgraph)$type <- "desired"
  }
  else {
    pos.subgraph <- dnet::dNetInduce(ig, pos.nodes, knn = 0, remove.loops = T,
                               largest.comp = F)
    conn.comp.graph <- igraph::decompose.graph(pos.subgraph)
    score.comp <- unlist(lapply(lapply(conn.comp.graph,
                                       function(x) as.numeric(V(x)$score)), sum))
    ind_order <- order(score.comp, decreasing = T)
    conn.comp.graph <- conn.comp.graph[ind_order]
    score.comp <- score.comp[ind_order]
    for (i in 1:length(conn.comp.graph)) {
      conn.comp.graph[[i]]$score <- score.comp[i]
    }
    v.id <- seq(1, vcount(ig))
    names(v.id) <- V(ig)$name
    edgelist <- igraph::get.edgelist(ig, names = F)
    edgelist1 <- edgelist[, 1]
    edgelist2 <- edgelist[, 2]
    #===#
    ig.size <- vcount(ig)
    for (i in 1:length(conn.comp.graph)) {
      new.id <- ig.size + i
      for (j in v.id[V(conn.comp.graph[[i]])$name]) {
        edgelist1[edgelist1 == j] <- new.id
        edgelist2[edgelist2 == j] <- new.id
      }
    }
    #===#
    new.ids <- seq(vcount(ig) + 1, vcount(ig) + length(conn.comp.graph))
    new.names <- paste("cluster", seq(1:length(conn.comp.graph)),
                       sep = "")
    names(new.ids) <- new.names
    v.id <- c(v.id, new.ids)
    v.name <- names(v.id)
    names(v.name) <- v.id
    new.edgelist <- cbind(v.name[as.character(edgelist1)],
                          v.name[as.character(edgelist2)])
    mig <- igraph::graph.edgelist(new.edgelist, directed = F)
    mig <- igraph::simplify(mig, remove.loops = T, remove.multiple = T)
    node.score <- scores[V(mig)$name]
    names(node.score) <- V(mig)$name
    node.score.cluster <- sapply(conn.comp.graph, igraph::get.graph.attribute,
                                 "score")
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
    if (!igraph::is.connected(mig)) {
      decomp.graphs <- igraph::decompose.graph(mig)
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
    all.ids <- c()
    if (length(mst.cluster.id) == 1) {
      neg.node.ids.2 = c()
    }
    else {
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
      sub.mig <- dnet::dNetInduce(mig, V(mst)$name, knn = 0, remove.loops = T, largest.comp = F)
      neg.node.ids <- which(V(sub.mig)$score < 0)
      for (i in neg.node.ids) {
        tmp_nei <- igraph::neighbors(sub.mig, v = i)
        tmp_nei_meta <- grep("cluster", V(sub.mig)[tmp_nei]$name)
        V(sub.mig)[i]$clusters <- list(tmp_nei[tmp_nei_meta])
      }
      score.neg.nodes <- c()
      for (i in neg.node.ids) {
        if (!is.na(V(sub.mig)[i]$clusters[1])) {
          borders <- c(i, V(sub.mig)[i]$clusters)
          borders <- unlist(borders)
          score.neg.nodes <- c(score.neg.nodes, sum(V(sub.mig)[borders]$score))
        }
        else {
          score.neg.nodes <- c(score.neg.nodes, V(sub.mig)[i]$score)
        }
      }
      neg.node.ids.2 <- neg.node.ids[score.neg.nodes > 0]
    }
    if (length(neg.node.ids.2) == 0) {
      tmp <- unlist(strsplit(names(node.score.cluster)[which.max(node.score.cluster)], "cluster"))
      ttmp <- as.numeric(matrix(tmp, nrow = 2)[2, ])
      tmp_nodes <- unlist(lapply(conn.comp.graph, igraph::get.vertex.attribute, "name")[ttmp])
      subgraph <- dnet::dNetInduce(ig, tmp_nodes, knn = 0, remove.loops = F,
                             largest.comp = T)
      V(subgraph)$type <- "desired"
    }
    else {
      subg <- dnet::dNetInduce(sub.mig, V(sub.mig)[neg.node.ids.2]$name,
                         knn = 0, remove.loops = T, largest.comp = T)
      mst.subg <- igraph::minimum.spanning.tree(subg, E(subg)$weight)
      getPathScore <- function(path, graph1, graph2) {
        s1 <- V(graph1)[path]$score
        tmp <- unique(unlist(V(graph1)[path]$clusters))
        s2 <- V(graph2)[tmp]$score
        sum(c(s1, s2))
      }
      max.score <- 0
      best.path <- c()
      for (i in 1:vcount(mst.subg)) {
        path <- igraph::all_shortest_paths(mst.subg, from = V(mst.subg)[i])
        path.score <- unlist(lapply(path$res, getPathScore,
                                    graph1 = mst.subg, graph2 = sub.mig))
        best.pos <- which.max(path.score)
        if (path.score[[best.pos]] > max.score) {
          best.path <- path$res[[best.pos]]
          max.score <- path.score[[best.pos]]
        }
      }
      if (length(best.path) != 1) {
        cluster.list <- V(mst.subg)[best.path]$clusters
        names.list <- as.character(1:length(cluster.list))
        names(cluster.list) <- names.list
        names(best.path) <- names.list
        for (i in names.list) {
          res <- lapply(cluster.list, intersect, cluster.list[[i]])
          if (length(intersect(unlist(cluster.list[as.character(which(as.numeric(names.list) < as.numeric(i)))]),
                               unlist(cluster.list[as.character(which(as.numeric(names.list) > as.numeric(i)))]))) > 0) {
            if (length(setdiff(res[[i]], unique(unlist(res[names(res) != i])))) == 0) {
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
      tmp_meta_nodes <- unlist(lapply(conn.comp.graph,
                                      igraph::get.vertex.attribute, "name")[ttmp])
      tmp_border_nodes <- V(mst.subg)[best.path]$name
      tmp_nodes <- c(tmp_border_nodes, tmp_meta_nodes)
      subgraph <- dnet::dNetInduce(ig, tmp_nodes, knn = 0, remove.loops = F,
                             largest.comp = T)
      type <- rep("desired", vcount(subgraph))
      names(type) <- V(subgraph)$name
      type[tmp_border_nodes[tmp_border_nodes %in% names(type)]] <- "linker"
      V(subgraph)$type <- type
    }
  }
  return(subgraph)
}
# EDIT 1/11: added examples, added description, added details


#' @title Exponential decay schedule for filtering rate
#'
#' @description Exponential decay schedule for filtering rate. Used to shift the raw RWR scores
#'
#' @param eta0 (numeric): starting filtering rate
#' @param d (numeric): decay parameter
#'
#' @return vector of filtering rates for each iteration of AMEND
#'
#' @examples
#' rates = AMEND:::exp.filtering.rate(eta0 = 0.5, d = 0.2)
#' plot(1:100, rates)
#'
exp.filtering.rate <- function(eta0 = 0.5, d = 0.1){
  iterations <- 1:100
  return( eta0*exp(-d*(iterations - 1)) )
}
# EDIT 1/11: added examples, added description


#' @title Sets the decay value
#'
#' @description Finds the max decay value for an exponential filtering rate schedule that will
#' allow the algorithm to arrive at a graph of size n.
#'
#' @param N (integer): Size of graph
#' @param eta0 (numeric): starting filtering rate
#' @param n (integer): approximate size of final module
#'
#' @return decay value, numeric
#'
#' @examples
#' d = AMEND:::find.decay(N = 1000, eta0 = 0.5, n = 50)
#' d
#'
find.decay <- function(N, eta0 = 0.5, n = 25){
  iter <- 1:100
  d.grid <- seq(0.01, 1, 0.01)
  decay <- c()

  for(j in 1:length(d.grid)){
    under <- 0
    perc <- exp.filtering.rate(eta0 = eta0, d = d.grid[j])

    d1 = N
    for(i in iter){
      d2 = round(d1 * (1 - perc[i] - 0.06), 0)
      if(d2 <= n){
        under <- 1
        break
      }
      if(d1 - d2 <= 2){
        break
      }
      d1 = d2
    }
    if(under){
      decay <- c(decay, d.grid[j])
    }else break
  }

  if(length(decay) == 0){
    return(NA)
  }else return(max(decay))
}
# EDIT 1/11: added examples


#' @title Calculate the core-clustering coefficients of a graph
#'
#' @description `core_cc()` calculates the core-clustering coefficients of all nodes in a graph. This network concept was introduced by Bader and Hogue (2003)
#'
#' @details The core-clustering coefficient of a node is the density of the maximum k-core of the immediate neighborhood of that node. The k-core of a graph is the maximal subgraph in which every node has degree >= k.
#'
#' @param g (igraph object): input graph. Must be connected
#'
#' @return Vector of core-clustering coefficients for each node
#'
#' @examples
#' core_clust_coefs = AMEND:::core_cc(glut4_graph)
#' head(core_clust_coefs)
#'
core_cc = function(g){
  if(!igraph::is_connected(g)) stop("Input graph must be connected")
  w = numeric(vcount(g))
  for(i in 1:vcount(g)){
    nei = as.numeric(igraph::neighbors(g, i))
    nborhood = igraph::induced_subgraph(g, c(nei, i))
    cores = igraph::coreness(nborhood)
    max_core = igraph::induced_subgraph(nborhood, which(cores == max(cores)))
    w[i] = igraph::edge_density(max_core)
  }
  return(w)
}
# EDIT 1/11: added examples, added description, added details

#=============#
# Testing ----
#=============#
# e = 0.6
# nc = 15
# eci.dir = "negative"
# og_adjM = readRDS("/Users/samboyd/Documents/GRA/Network Analysis/AMEND/Subnetwork Data/glut4_adjM.RDS")
# og_node_scores = readRDS("/Users/samboyd/Documents/GRA/Network Analysis/AMEND/Subnetwork Data/glut4_node_scores.RDS")
#
# subnet = amend(eta = e, adj_matrix = og_adjM, node_scores = og_node_scores, eci.direction = eci.dir, n = nc)


