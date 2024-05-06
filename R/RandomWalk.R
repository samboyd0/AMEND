#' @importFrom igraph vcount V V<- E E<- vertex_attr vertex_attr<-
#' @importFrom foreach %dopar%

#' @title Create a transition matrix from the adjacency matrix of a graph
#'
#' @description
#' `transition_matrix()` creates a transition matrix from the adjacency matrix of a graph for input into random walk with restart (RWR).
#'
#' @details
#' There are two normalization schemes available.
#' * degree: Divides each column of the adjacency matrix by the column total
#' * modified_degree: First creates a modified adjacency matrix by left and right multiplying the adjacency matrix by a diagonal matrix raised to the power of -k, whose diagonal elements are the column sums of the original adjacency matrix. Then proceed to column-normalize this modified adjacency matrix.
#'
#' For multiplex/heterogeneous networks, the normalization process requires more steps. First, all multiplex and layer-to-layer bipartite matrices are column-normalized independently. For the multiplex component, this involves normalizing each inter- and intra-layer matrix independently such that their columns sum to one. Then _switch.layer.prob_ is applied to the columns of the entire multiplex adjacency matrix such that these now sum to one. Finally, _jump.prob_ is applied to the columns of the supra-adjacency matrix (i.e., matrix of the whole integrated graph) such that these sum to one.
#'
#' @param adjM Adjacency matrix
#' @param norm Normalization method
#' @param k Value between 0 and 1. When norm = "modified_degree", the adjacency matrix is first left and right multiplied by a diagonal matrix of node degrees raised to the power -k. As k increases, edge weights are penalized more for the degrees of their adjacent nodes.
#' @param heterogeneous Logical. If TRUE, graph is considered heterogeneous (more than one distinct node type, e.g., proteins and metabolites)
#' @param multiplex Logical. If true, graph is assumed to contain multiplex components.
#' @param jump.prob A named vector, or NULL. Probability of random walker jumping from one component of graph to another in RWR. Only used when heterogeneous=TRUE.
#' @param switch.layer.prob A named list of named vectors, or NULL. Probability of random walker to switch from current layer in a multiplex to another layer in same component. List element names correspond to multiplex components, and vector names correspond to layers within a multiplex.
#' @param brw.attr A numeric vector or NULL. Biased random walk vertex attribute values. Should be non-negative, with values greater (lesser) than 1 increasing (decreasing) transition probabilities to a node in RWR. If NULL, all nodes are given a value of 1.
#' @param degree.bias List or NULL. NULL (default) for no degree bias adjustment. List names should include 'component', which is a character string of length greater than or equal to 1 specifying the components or layers to which the degree bias adjustment will be applied, and 'method', which is a character string of 'BS' (bistochastic scaling), 'IN' (inflation-normalization), or 'SDS' (stationary distribution scaling).
#' @param in.parallel Logical. Run computations in parallel.
#' @param n.cores number of cores to use in parallel computations
#'
#' @return transition matrix
#'
#' @seealso [run_AMEND()], [RWR()], [create_integrated_graph()]
#'
#' @examples
#'
#' library(igraph)
#'
#' # Creating an adjacency matrix and graph
#' adjm = Matrix::Matrix(c(0, 1, 0, 0, 0, 0, 0, 0,
#'                       1, 0, 1, 1, 1, 0, 0, 0,
#'                       0, 1, 0, 0, 1, 1, 0, 0,
#'                       0, 1, 0, 0, 0, 0, 1, 0,
#'                       0, 1, 1, 0, 0, 0, 1, 1,
#'                       0, 0, 1, 0, 0, 0, 0, 0,
#'                       0, 0, 0, 1, 1, 0, 0, 0,
#'                       0, 0, 0, 0, 1, 0, 0, 0), nrow = 8)
#' g = graph_from_adjacency_matrix(adjm, mode = 'undirected')
#'
#' adj_norm = transition_matrix(adjM = adjm, norm = "degree")
#'
#' @export
transition_matrix <- function(adjM, norm = c("degree", "modified_degree"), k = 0.5, heterogeneous = FALSE, multiplex = FALSE,
                              jump.prob, switch.layer.prob, brw.attr = NULL, degree.bias = NULL, in.parallel = FALSE, n.cores){
  # Extract dimnames
  adj.dimnames = dimnames(adjM)

  # Coerce to a general, column-compressed sparse matrix from Matrix package
  if(!"dgCMatrix" %in% class(adjM)) adjM = methods::as(methods::as(methods::as(adjM, "dMatrix"), "generalMatrix"), "CsparseMatrix")

  # Apply Biased Random Walk values to adjacency matrix before normalization
  if(is.null(brw.attr)) brw.attr = rep(1, nrow(adjM))
  adjM = Matrix::Diagonal(x = brw.attr) %*% adjM
  dimnames(adjM) = adj.dimnames

  # Extract node type from rownames
  if(multiplex || heterogeneous) node_type = get.type(rownames(adjM),2)

  if(!heterogeneous && !multiplex){
    nadjM = column_normalize(adjM, norm, k, dimnames(adjM))
    if(!is.null(degree.bias)){
      if(degree.bias$method == 'BS') nadjM = bistochastic_scaling(nadjM)
      if(degree.bias$method == 'IN') nadjM = inflate_normalize(nadjM)
    }
  }else if(heterogeneous && !multiplex){
    # Normalize monoplex components
    uniq.types = unique(node_type)
    monoplex.comps = uniq.types[!grepl("_", uniq.types)]
    mono.adj = vector("list", length(monoplex.comps)); names(mono.adj) = monoplex.comps
    for(i in seq_along(mono.adj)){
      id = which(get.type(rownames(adjM), 3) %in% monoplex.comps[i])
      mono.adj[[i]] = adjM[id,id,drop=FALSE]
    }

    if(!in.parallel){
      mono.norm = vector("list", length(monoplex.comps)); names(mono.norm) = monoplex.comps
      for(i in seq_along(mono.adj)){
        new.mat = column_normalize(mono.adj[[i]], norm, k, dimnames(mono.adj[[i]]))
        mono.norm[[i]] = sum2one(new.mat)
        if(!is.null(degree.bias) && names(mono.norm)[i] %in% degree.bias$component){
          if(degree.bias$method == 'BS') mono.norm[[i]] = bistochastic_scaling(mono.norm[[i]])
          if(degree.bias$method == 'IN') mono.norm[[i]] = inflate_normalize(mono.norm[[i]])
        }
      }
    }
    if(in.parallel){
      cl = parallel::makeForkCluster(n.cores, outfile = "")
      doParallel::registerDoParallel(cl)
      mono.norm = foreach::foreach(i = 1:length(monoplex.comps), .verbose = FALSE, .packages = c("igraph", "Matrix"), .export = NULL, .noexport = NULL) %dopar% {
        new.mat = column_normalize(mono.adj[[i]], norm, k, dimnames(mono.adj[[i]]))
        new.mat = sum2one(new.mat)
        if(!is.null(degree.bias) && monoplex.comps[i] %in% degree.bias$component){
          if(degree.bias$method == 'BS') new.mat = bistochastic_scaling(new.mat)
          if(degree.bias$method == 'IN') new.mat = inflate_normalize(new.mat)
        }
        new.mat
      }
      parallel::stopCluster(cl)
      names(mono.norm) = monoplex.comps
    }

    # Normalize off-diagonal sub-matrices (and sub-sub matrices), then combine adjacency matrices
    uniq.types = unique(extract_string(node_type, "_", pos=1))
    new.mat = mono.norm[[1]]
    for(i in 2:length(uniq.types)){ # Outside (row)
      for(j in 1:(i-1)){ # Inside (col)
        nt.r = get.type(rownames(mono.norm[[i]]),2)
        uniq.nt.r = unique(nt.r)
        nt.c = get.type(colnames(mono.norm[[j]]),2)
        uniq.nt.c = unique(nt.c)
        bottom.left.tmp = Matrix::Matrix(data = 0, nrow = nrow(mono.norm[[i]]), ncol = ncol(mono.norm[[j]]), dimnames = list(rownames(mono.norm[[i]]), colnames(mono.norm[[j]])), sparse = TRUE)
        for(ii in seq_along(uniq.nt.r)){ # Row-wise
          for(jj in seq_along(uniq.nt.c)){ # Col-wise
            tmp = adjM[which(get.type(rownames(adjM),2) == uniq.nt.r[ii]), which(get.type(colnames(adjM),2) == uniq.nt.c[jj]),drop=FALSE]
            tmp.dimnames = list(dimnames(adjM)[[1]][which(get.type(rownames(adjM),2) == uniq.nt.r[ii])], dimnames(adjM)[[2]][which(get.type(colnames(adjM),2) == uniq.nt.c[jj])])
            tmp = column_normalize(tmp, norm, k, tmp.dimnames)
            bottom.left.tmp[which(nt.r == uniq.nt.r[ii]), which(nt.c == uniq.nt.c[jj])] = tmp
          }
        }
        nt.r = get.type(rownames(mono.norm[[j]]),2)
        nt.c = get.type(colnames(mono.norm[[i]]),2)
        uniq.nt.r = unique(nt.r)
        uniq.nt.c = unique(nt.c)
        top.right.tmp = Matrix::Matrix(data = 0, nrow = nrow(mono.norm[[j]]), ncol = ncol(mono.norm[[i]]), dimnames = list(rownames(mono.norm[[j]]), colnames(mono.norm[[i]])), sparse = TRUE)
        for(ii in seq_along(uniq.nt.r)){ # Row-wise
          for(jj in seq_along(uniq.nt.c)){ # Col-wise
            tmp = adjM[which(get.type(rownames(adjM),2) == uniq.nt.r[ii]), which(get.type(colnames(adjM),2) == uniq.nt.c[jj]),drop=FALSE]
            tmp.dimnames = list(dimnames(adjM)[[1]][which(get.type(rownames(adjM),2) == uniq.nt.r[ii])], dimnames(adjM)[[2]][which(get.type(colnames(adjM),2) == uniq.nt.c[jj])])
            tmp = column_normalize(tmp, norm, k, tmp.dimnames)
            top.right.tmp[which(nt.r == uniq.nt.r[ii]), which(nt.c == uniq.nt.c[jj])] = tmp
          }
        }
        if(j == 1){
          top.right.mat = top.right.tmp
          bottom.left.mat = bottom.left.tmp
        }else{
          top.right.mat = rbind(top.right.mat, top.right.tmp)
          bottom.left.mat = cbind(bottom.left.mat, bottom.left.tmp)
        }
      }
      new.mat = cbind(rbind(new.mat, bottom.left.mat), rbind(top.right.mat, mono.norm[[i]]))
    }

    ## Apply jump.prob (which will be a named vector)
    ntl = get.type(rownames(new.mat),2) # node type w/ layers
    nt = get.type(rownames(new.mat),3) # node types
    for(j in seq_len(ncol(new.mat))){
      x.ids = (new.mat@p[j] + 1):new.mat@p[j+1]
      r.ids = new.mat@i[x.ids] + 1 # row ids of non-zero elements in column j
      nnt = extract_string(unique(ntl[r.ids][nt[r.ids] != nt[j]]), "_", pos=1)
      nL = table(nnt) # number of layers
      ntc = length(nL) # Number of other node type connections for node j
      current = nt[j] # current node type
      other = nt[r.ids] # other node types
      tmp = ifelse(all(other != current), 1, jump.prob[current])
      if(ntc == 0) jp = 1 else jp = ifelse(other == current, 1 - jump.prob[current], tmp / ( ntc * nL[match(other, names(nL))] ) )
      new.mat@x[x.ids] = new.mat@x[x.ids] * jp
    }
    nadjM = new.mat
  }else if(multiplex && !heterogeneous){
    # Normalize multiplex components
    ## Get Multiplex components as a list of adjacency matrices
    uniq.types = unique(node_type)
    multiplex.comps = unique(extract_string(uniq.types[grepl("_", uniq.types)], "_", pos=1)) # Should be length 1
    id = which(get.type(rownames(adjM), 3) %in% multiplex.comps)
    multi.adj = adjM[id,id,drop=FALSE]

    ## Number of layers
    L = length(unique(get.type(rownames(multi.adj), 2)))
    ## IDs of columns for each layer as a list
    uniq.l = unique(get.type(rownames(multi.adj),4))
    l.ids = vector("list", L)
    for(j in seq_along(l.ids)){
      l.ids[[j]] = which(get.type(rownames(multi.adj),4) == uniq.l[j])
    }
    ## Number of layers each node connects to
    nL = numeric(nrow(multi.adj)); names(nL) = rownames(multi.adj)
    for(j in seq_len(ncol(multi.adj))){
      id = multi.adj@i[(multi.adj@p[j] + 1):multi.adj@p[j+1]] + 1 # row ids of non-zero elements in column j
      nnt = get.type(rownames(multi.adj)[id],2)
      nnt = nnt[nnt != get.type(rownames(multi.adj)[j],2)]
      nL[j] = length(unique(nnt))
    }
    ## Normalize each sub-matrix of each multiplex component
    new.mat = Matrix::Matrix(data = 0, nrow = nrow(multi.adj), ncol = ncol(multi.adj), dimnames = dimnames(multi.adj), sparse = TRUE)
    for(j in seq_len(L)){ # Row-wise
      for(ii in seq_len(L)){ # Col-wise
        tmp = multi.adj[l.ids[[j]], l.ids[[ii]],drop=FALSE]
        tmp.dimnames = list(dimnames(multi.adj)[[1]][l.ids[[j]]], dimnames(multi.adj)[[2]][l.ids[[ii]]])
        tmp = column_normalize(tmp, norm, k, tmp.dimnames)
        if(j == ii & !is.null(degree.bias) && paste(multiplex.comps, uniq.l[j], sep='_') %in% degree.bias$component){
          if(degree.bias$method == 'BS') tmp = bistochastic_scaling(tmp)
          if(degree.bias$method == 'IN') tmp = inflate_normalize(tmp)
        }
        new.mat[l.ids[[j]], l.ids[[ii]]] = tmp
      }
    }
    ## Apply switch.layer.prob (which will be a named list of named vectors)
    nt = get.type(rownames(new.mat),2)
    slp = switch.layer.prob[[multiplex.comps]]
    for(j in seq_len(ncol(new.mat))){
      x.ids = (new.mat@p[j] + 1):new.mat@p[j+1]
      current.layer = nt[j]
      other.layers = nt[new.mat@i[x.ids] + 1]
      tmp = ifelse(all(other.layers != current.layer), 1, slp[current.layer])
      if(nL[j] == 0) switch.probs = 1 else switch.probs = ifelse(other.layers == current.layer, 1 - slp[current.layer], tmp / nL[j])
      new.mat@x[x.ids] = new.mat@x[x.ids] * switch.probs
    }
    nadjM = sum2one(new.mat)
  }else if(multiplex && heterogeneous){
    # Normalize multiplex components
    ## Get Multiplex components as a list of adjacency matrices
    uniq.types = unique(node_type)
    multiplex.comps = unique(extract_string(uniq.types[grepl("_", uniq.types)], "_", pos=1))
    multi.adj = vector("list", length(multiplex.comps)); names(multi.adj) = multiplex.comps
    for(i in seq_along(multi.adj)){
      id = which(get.type(rownames(adjM), 3) %in% multiplex.comps[i])
      multi.adj[[i]] = adjM[id,id,drop=FALSE]
    }

    if(!in.parallel){
      multi.norm = vector("list", length(multiplex.comps)); names(multi.norm) = multiplex.comps
      for(i in seq_along(multi.adj)){ # For each multiplex component
        # Column sums to determine columns of all zeros
        cs = Matrix::colSums(multi.adj[[i]])
        ## Number of layers
        L = length(unique(get.type(rownames(multi.adj[[i]]), 2)))
        ## IDs of columns for each layer as a list
        uniq.l = unique(get.type(rownames(multi.adj[[i]]),4))
        l.ids = vector("list", L)
        for(j in seq_along(l.ids)){
          l.ids[[j]] = which(get.type(rownames(multi.adj[[i]]),4) == uniq.l[j])
        }
        ## Number of other layers each node connects to
        nL = numeric(nrow(multi.adj[[i]])); names(nL) = rownames(multi.adj[[i]])
        for(j in which(cs != 0)){ # seq_len(ncol(multi.adj[[i]]))
          id = multi.adj[[i]]@i[(multi.adj[[i]]@p[j] + 1):multi.adj[[i]]@p[j+1]] + 1 # row ids of non-zero elements in column j
          nnt = get.type(rownames(multi.adj[[i]])[id],2)
          nnt = nnt[nnt != get.type(rownames(multi.adj[[i]])[j],2)]
          nL[j] = length(unique(nnt))
        }
        ## Normalize each sub-matrix of each multiplex component
        new.mat = Matrix::Matrix(data = 0, nrow = nrow(multi.adj[[i]]), ncol = ncol(multi.adj[[i]]), dimnames = dimnames(multi.adj[[i]]), sparse = TRUE)
        for(j in seq_len(L)){ # Row-wise
          for(ii in seq_len(L)){ # Col-wise
            tmp = multi.adj[[i]][l.ids[[j]], l.ids[[ii]],drop=FALSE]
            tmp.dimnames = list(dimnames(multi.adj[[i]])[[1]][l.ids[[j]]], dimnames(multi.adj[[i]])[[2]][l.ids[[ii]]])
            tmp = column_normalize(tmp, norm, k, tmp.dimnames)
            if(j == ii & !is.null(degree.bias) && paste(multiplex.comps[i], uniq.l[j], sep='_') %in% degree.bias$component){
              if(degree.bias$method == 'BS') tmp = bistochastic_scaling(tmp)
              if(degree.bias$method == 'IN') tmp = inflate_normalize(tmp)
            }
            new.mat[l.ids[[j]], l.ids[[ii]]] = tmp
          }
        }
        ## Apply switch.layer.prob (which will be a named list of named vectors)
        nt = get.type(rownames(new.mat),2)
        slp = switch.layer.prob[[names(multi.adj)[i]]]
        for(j in which(cs != 0)){
          x.ids = (new.mat@p[j] + 1):new.mat@p[j+1]
          current.layer = nt[j]
          other.layers = nt[new.mat@i[x.ids] + 1]
          tmp = ifelse(all(other.layers != current.layer), 1, slp[current.layer])
          if(nL[j] == 0) switch.probs = 1 else switch.probs = ifelse(other.layers == current.layer, 1 - slp[current.layer], tmp / nL[j])
          new.mat@x[x.ids] = new.mat@x[x.ids] * switch.probs
        }
        multi.norm[[i]] = sum2one(new.mat)
      }
    }
    if(in.parallel){
      cl = parallel::makeForkCluster(n.cores, outfile = "")
      doParallel::registerDoParallel(cl)
      multi.norm = foreach::foreach(i = 1:length(multi.adj), .verbose = FALSE, .packages = c("igraph", "Matrix"), .export = NULL, .noexport = NULL) %dopar% {
        # Column sums to determine columns of all zeros
        cs = Matrix::colSums(multi.adj[[i]])
        ## Number of layers
        L = length(unique(get.type(rownames(multi.adj[[i]]), 2)))
        ## IDs of columns for each layer as a list
        uniq.l = unique(get.type(rownames(multi.adj[[i]]),4))
        l.ids = vector("list", L)
        for(j in seq_along(l.ids)){
          l.ids[[j]] = which(get.type(rownames(multi.adj[[i]]),4) == uniq.l[j])
        }
        ## Number of other layers each node connects to
        nL = numeric(nrow(multi.adj[[i]])); names(nL) = rownames(multi.adj[[i]])
        for(j in which(cs != 0)){
          id = multi.adj[[i]]@i[(multi.adj[[i]]@p[j] + 1):multi.adj[[i]]@p[j+1]] + 1 # row ids of non-zero elements in column j
          nnt = get.type(rownames(multi.adj[[i]])[id],2)
          nnt = nnt[nnt != get.type(rownames(multi.adj[[i]])[j],2)]
          nL[j] = length(unique(nnt))
        }
        ## Normalize each sub-matrix of each multiplex component
        new.mat = Matrix::Matrix(data = 0, nrow = nrow(multi.adj[[i]]), ncol = ncol(multi.adj[[i]]), dimnames = dimnames(multi.adj[[i]]), sparse = TRUE)
        for(j in seq_len(L)){ # Row-wise
          for(ii in seq_len(L)){ # Col-wise
            tmp = multi.adj[[i]][l.ids[[j]], l.ids[[ii]],drop=FALSE]
            tmp.dimnames = list(dimnames(multi.adj[[i]])[[1]][l.ids[[j]]], dimnames(multi.adj[[i]])[[2]][l.ids[[ii]]])
            tmp = column_normalize(tmp, norm, k, tmp.dimnames)
            if(j == ii & !is.null(degree.bias) && paste(multiplex.comps[i], uniq.l[j], sep='_') %in% degree.bias$component){
              if(degree.bias$method == 'BS') tmp = bistochastic_scaling(tmp)
              if(degree.bias$method == 'IN') tmp = inflate_normalize(tmp)
            }
            new.mat[l.ids[[j]], l.ids[[ii]]] = tmp
          }
        }
        ## Apply switch.layer.prob (which will be a named list of named vectors)
        nt = get.type(rownames(new.mat),2)
        slp = switch.layer.prob[[names(multi.adj)[i]]]
        for(j in which(cs != 0)){
          x.ids = (new.mat@p[j] + 1):new.mat@p[j+1]
          current.layer = nt[j]
          other.layers = nt[new.mat@i[x.ids] + 1]
          tmp = ifelse(all(other.layers != current.layer), 1, slp[current.layer])
          if(nL[j] == 0) switch.probs = 1 else switch.probs = ifelse(other.layers == current.layer, 1 - slp[current.layer], tmp / nL[j])
          new.mat@x[x.ids] = new.mat@x[x.ids] * switch.probs
        }
        sum2one(new.mat)
      }
      parallel::stopCluster(cl)
      names(multi.norm) = multiplex.comps
    }

    # Normalize monoplex components
    if(any(!grepl("_", uniq.types))){
      monoplex.comps = uniq.types[!grepl("_", uniq.types)]
      mono.adj = vector("list", length(monoplex.comps)); names(mono.adj) = monoplex.comps
      for(i in seq_along(mono.adj)){
        id = which(get.type(rownames(adjM), 3) %in% monoplex.comps[i])
        mono.adj[[i]] = adjM[id,id,drop=FALSE]
      }
      if(!in.parallel){
        mono.norm = vector("list", length(monoplex.comps)); names(mono.norm) = monoplex.comps
        for(i in seq_along(mono.adj)){
          new.mat = column_normalize(mono.adj[[i]], norm, k, dimnames(mono.adj[[i]]))
          mono.norm[[i]] = sum2one(new.mat)
          if(!is.null(degree.bias) && names(mono.norm)[i] %in% degree.bias$component){
            if(degree.bias$method == 'BS') mono.norm[[i]] = bistochastic_scaling(mono.norm[[i]])
            if(degree.bias$method == 'IN') mono.norm[[i]] = inflate_normalize(mono.norm[[i]])
          }
        }
      }
      if(in.parallel){
        cl = parallel::makeForkCluster(n.cores, outfile = "")
        doParallel::registerDoParallel(cl)
        mono.norm = foreach::foreach(i = 1:length(monoplex.comps), .verbose = FALSE, .packages = c("igraph", "Matrix"), .export = NULL, .noexport = NULL) %dopar% {
          new.mat = column_normalize(mono.adj[[i]], norm, k, dimnames(mono.adj[[i]]))
          new.mat = sum2one(new.mat)
          if(!is.null(degree.bias) && monoplex.comps[i] %in% degree.bias$component){
            if(degree.bias$method == 'BS') new.mat = bistochastic_scaling(new.mat)
            if(degree.bias$method == 'IN') new.mat = inflate_normalize(new.mat)
          }
          new.mat
        }
        parallel::stopCluster(cl)
        names(mono.norm) = monoplex.comps
      }
      adj.comps = c(multi.norm, mono.norm)
    }else adj.comps = multi.norm

    # Normalize off-diagonal sub-matrices (and sub-sub matrices), then combine adjacency matrices
    uniq.types = unique(extract_string(node_type, "_", pos=1))
    new.mat = adj.comps[[1]]
    for(i in 2:length(uniq.types)){ # Outside (row)
      for(j in 1:(i-1)){ # Inside (col)
        nt.r = get.type(rownames(adj.comps[[i]]),2) ### i subscript out of bounds! (iteration 2 of run_AMEND)
        uniq.nt.r = unique(nt.r)
        nt.c = get.type(colnames(adj.comps[[j]]),2)
        uniq.nt.c = unique(nt.c)
        bottom.left.tmp = Matrix::Matrix(data = 0, nrow = nrow(adj.comps[[i]]), ncol = ncol(adj.comps[[j]]), dimnames = list(rownames(adj.comps[[i]]), colnames(adj.comps[[j]])), sparse = TRUE)
        for(ii in seq_along(uniq.nt.r)){ # Row-wise
          for(jj in seq_along(uniq.nt.c)){ # Col-wise
            tmp = adjM[which(get.type(rownames(adjM),2) == uniq.nt.r[ii]), which(get.type(colnames(adjM),2) == uniq.nt.c[jj]),drop=FALSE]
            tmp.dimnames = list(dimnames(adjM)[[1]][which(get.type(rownames(adjM),2) == uniq.nt.r[ii])], dimnames(adjM)[[2]][which(get.type(colnames(adjM),2) == uniq.nt.c[jj])])
            tmp = column_normalize(tmp, norm, k, tmp.dimnames)
            bottom.left.tmp[which(nt.r == uniq.nt.r[ii]), which(nt.c == uniq.nt.c[jj])] = tmp
          }
        }
        nt.r = get.type(rownames(adj.comps[[j]]),2)
        nt.c = get.type(colnames(adj.comps[[i]]),2)
        uniq.nt.r = unique(nt.r)
        uniq.nt.c = unique(nt.c)
        top.right.tmp = Matrix::Matrix(data = 0, nrow = nrow(adj.comps[[j]]), ncol = ncol(adj.comps[[i]]), dimnames = list(rownames(adj.comps[[j]]), colnames(adj.comps[[i]])), sparse = TRUE)
        for(ii in seq_along(uniq.nt.r)){ # Row-wise
          for(jj in seq_along(uniq.nt.c)){ # Col-wise
            tmp = adjM[which(get.type(rownames(adjM),2) == uniq.nt.r[ii]), which(get.type(colnames(adjM),2) == uniq.nt.c[jj]),drop=FALSE]
            tmp.dimnames = list(dimnames(adjM)[[1]][which(get.type(rownames(adjM),2) == uniq.nt.r[ii])], dimnames(adjM)[[2]][which(get.type(colnames(adjM),2) == uniq.nt.c[jj])])
            tmp = column_normalize(tmp, norm, k, tmp.dimnames)
            top.right.tmp[which(nt.r == uniq.nt.r[ii]), which(nt.c == uniq.nt.c[jj])] = tmp
          }
        }
        if(j == 1){
          top.right.mat = top.right.tmp
          bottom.left.mat = bottom.left.tmp
        }else{
          top.right.mat = rbind(top.right.mat, top.right.tmp)
          bottom.left.mat = cbind(bottom.left.mat, bottom.left.tmp)
        }
      }
      new.mat = cbind(rbind(new.mat, bottom.left.mat), rbind(top.right.mat, adj.comps[[i]]))
    }

    ## Apply jump.prob (which will be a named vector)
    ntl = get.type(rownames(new.mat),2) # node types w/ layers (e.g., 'prot_1', 'prot_2', 'gene', ...)
    nt = get.type(rownames(new.mat),3) # node types (e.g., 'prot', 'gene', ...)
    for(j in seq_len(ncol(new.mat))){
      x.ids = (new.mat@p[j] + 1):new.mat@p[j+1]
      r.ids = new.mat@i[x.ids] + 1 # row ids of non-zero elements in column j (neighbor's of j)
      ont = extract_string(unique(ntl[r.ids][nt[r.ids] != nt[j]]), "_", pos=1) # (e.g., 'prot', 'prot', 'gene', ...)
      nL = table(ont) # number of layers of each other component that node j is connected to
      ntc = length(nL) # Number of connections with other node types for node j
      current = nt[j] # current node type
      other = nt[r.ids] # other node types
      # If node j is only connected to nodes of another type, then jump prob = 1
      tmp.jp = ifelse(all(other != current), 1, jump.prob[current])
      if(ntc == 0) jp = 1 else jp = ifelse(other == current, 1 - jump.prob[current], tmp.jp / ( ntc * nL[match(other, names(nL))] ) ) # If ntc == 0, then all(other == current)
      new.mat@x[x.ids] = new.mat@x[x.ids] * jp
    }
    nadjM = new.mat
  }
  return(nadjM)
}

#' @title Random walk with restart (RWR) procedure
#'
#' @description
#' `RWR()` implements the RWR procedure. This simulates random walkers on a graph, with a certain probability of returning to the seeds nodes given by the restart parameter.
#'
#' @details
#' Seed values must be non-negative. For heterogeneous-multiplex graphs, the seed vector is a weighted concatenation of seed vectors from each component. Weights are given by _net.weight_ and must sum to one. These component-wise seed vectors are themselves weighted concatenations of seed vectors from each layer in the component, with each layer-wise seed vector independently normalized to sum to one. Weights are given by _layer.weight_ and must sum to one within each component. After the construction of the transition matrix and seed vector, RWR is implemented as usual through iterative matrix multiplication until convergence.
#'
#' @param nadjM Normalized adjacency matrix
#' @param setSeeds A named vector of seed values
#' @param restart Restart parameter
#' @param heterogeneous Logical. If TRUE, graph is considered heterogeneous (more than one distinct node type, e.g., proteins and metabolites)
#' @param multiplex Logical. If true, graph is assumed to contain multiplex components.
#' @param net.weight A named vector, or NULL. Relative weight given to nodes of a component of graph, applied to seed vector in RWR. Only used when heterogeneous=TRUE.
#' @param layer.weight A named list of named vectors, or NULL. Relative weight given to nodes of a layer of a component of graph, applied to seed vector in RWR. List element names correspond to multiplex components, and vector names correspond to layers within a multiplex.
#' @param stop_step Maximum number of iterations of RWR.
#'
#' @return matrix of diffusion scores from RWR
#'
#' @seealso [run_AMEND()], [transition_matrix()], [create_integrated_graph()]
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
#' adj_norm = Matrix::Matrix(adjm %*% diag(1 / degree(g)), dimnames = dimnames(adjm), sparse = TRUE)
#'
#' # Creating a named vector of seed values
#' seeds = runif(8)
#' names(seeds) = 1:8
#'
#' rwr = RWR(nadjM = adj_norm, setSeeds = seeds, restart = 0.8)
#'
#' @export
RWR <- function (nadjM, setSeeds = NULL, restart = 0.75, heterogeneous = FALSE, multiplex = FALSE, net.weight, layer.weight, stop_step = 50){
  # NB: In the context of run_AMEND(), nadjM and setSeeds will be in the same order
  if (is.null(restart) || is.na(restart) || restart < 0 || restart > 100) {
    r <- 0.75
  }else if (restart > 1 && restart < 100) {
    r <- restart/100
  }else {
    r <- restart
  }
  stop_delta <- 1e-06
  # stop_step <- 50
  if (is.null(setSeeds)){
    # stop("setSeeds must be non-NULL")
    setSeeds = matrix(1, ncol = 1, nrow = nrow(nadjM), dimnames = list(rownames(nadjM)))
  }
  node_type = get.type(rownames(nadjM),2)
  if (is.data.frame(setSeeds)) {
    data <- as.matrix(setSeeds)
  }else if (is.vector(setSeeds)) {
    data <- as.matrix(setSeeds, ncol = 1)
  }else data <- setSeeds
  if (is.null(rownames(data))) {
    stop("The function must require the row names of the input setSeeds.\n")
  }else if (any(is.na(rownames(data)))) {
    warning("setSeeds with NA as row names will be removed")
    data <- data[!is.na(rownames(data)), ,drop=FALSE]
    node_type = node_type[!is.na(rownames(data))]
  }
  if(multiplex && heterogeneous){
    layers = unique(node_type)
    p0 = numeric(length(node_type))
    for(i in seq_along(layers)){
      id = which(node_type == layers[i])
      if(grepl("_", layers[i])) tmp = layer.weight[[extract_string(layers[i], "_", pos=1)]][layers[i]] else tmp = 1
      p0[id] = sum2one(data[id, ,drop=FALSE]) * tmp
    }
    nt = unique(extract_string(layers, "_", pos=1))
    for(i in seq_along(nt)){
      id = which(extract_string(node_type, "_", pos=1) == nt[i])
      p0[id] = p0[id] * net.weight[nt[i]]
    }
    P0matrix = Matrix::Matrix(p0, ncol = 1, sparse = TRUE)
  }
  if(!multiplex && heterogeneous){
    nt = unique(node_type)
    p0 = numeric(length(node_type))
    for(i in seq_along(nt)){
      id = which(node_type == nt[i])
      p0[id] = sum2one(data[id, ,drop=FALSE]) * net.weight[nt[i]]
    }
    P0matrix = Matrix::Matrix(p0, ncol = 1, sparse = TRUE)
  }
  if(multiplex && !heterogeneous){
    layers = unique(node_type)
    p0 = numeric(length(node_type))
    for(i in seq_along(layers)){
      id = which(node_type == layers[i])
      p0[id] = sum2one(data[id, ,drop=FALSE]) * layer.weight[[extract_string(layers[i], "_", pos=1)]][layers[i]]
    }
    P0matrix = Matrix::Matrix(p0, ncol = 1, sparse = TRUE)
  }
  if(!multiplex && !heterogeneous){
    P0matrix <- sum2one(data)
  }
  if (restart == 1) {
    PTmatrix <- P0matrix
  }else {
    PTmatrix <- Matrix::Matrix(0, nrow = nrow(P0matrix), ncol = ncol(P0matrix), sparse = T)
    P0 <- P0matrix[, 1,drop=FALSE]
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
  rownames(PTmatrix) <- rownames(P0matrix)
  colnames(PTmatrix) <- colnames(P0matrix)
  invisible(PTmatrix)
}

#' @title Compute entropy of a discrete probability distribution
#'
#' @description
#' Computes entropy (with natural log) of a discrete probability distribution vector (i.e., sums to 1, values between 0 and 1).
#'
#' @param x Probability vector
#'
#' @returns Numeric value representing entropy
#'
entropy = function(x){
  tmp = ifelse(round(x, 10) == 0, 0, -x * log(x))
  sum(tmp)
}

#' @title Compute stationary distribution of a transition matrix
#'
#' @description
#' Compute the stationary distribution of a Markov chain represented by a transition matrix. This is the eigenvector associated with the absolute largest eigenvalue of the transition matrix (eigenvector normalized to sum to one).
#'
#' @param x Left-stochastic transition matrix
#'
#' @returns Probability vector of the stationary distribution.
#'
stationary.distr = function(x){
  # e = Re(RSpectra::eigs(A = x, k = 1, which = "LM")$vectors[,1])
  # tmp = e / sum(e)
  # ifelse(tmp < 0, 0, tmp)
  e = RWR(nadjM = x, restart = 0, stop_step = 100)
  e[,1]
}

#' @title Scale edge weights of adjacency matrix by degrees of adjacent nodes.
#'
#' @description
#' This functions creates a modified adjacency matrix. The adjacency matrix is left and right multiplied by two identical diagonal matrices whose diagonal elements are node degrees raised to the power -_k_.
#'
#' This decreases edge weights in proportion to the degrees of their adjacent nodes. In the context of column-normalization to obtain a transition matrix, this results in a biased random walk where the transition probability to a target node is inversely proportional to the target node's degree (all else being equal).
#'
#' This can aid in mitigating degree bias associated with PPI networks.
#'
#' @param adjM Adjacency matrix
#' @param k Exponent applied to node degrees. Should be positive, since degrees are raised to -_k_.
#'
#' @returns Modified adjacency matrix
#'
modify_adj_mat = function(adjM, k){
  # For symmetric matrices, inv.col.sum = inv.row.sum
  adj.dimnames = dimnames(adjM)
  inv.col.sum <- Matrix::Diagonal(x = Matrix::colSums(adjM)^(-k))
  inv.row.sum <- Matrix::Diagonal(x = Matrix::rowSums(adjM)^(-k))
  adjM.mod <- inv.row.sum %*% adjM %*% inv.col.sum
  dimnames(adjM.mod) = adj.dimnames
  return(adjM.mod)
}

#' @title Column-normalize an adjacency matrix
#'
#' @description
#' Given the normalization scheme, this computes a column-normalized adjacency matrix.
#'
#' When norm='degree', the columns of the adjacency matrix are divided by their column sums. When norm='modified_degree', the adjacency matrix is first modified to decrease edge weights in proportion to the degree of its adjacent nodes, then it is column-normalized as with norm='degree'.
#'
#' @param adjM Adjacency matrix
#' @param norm Normalization scheme. One of 'degree' or 'modified_degree'
#' @param k Exponent applied to node degrees. Should be positive, since degrees are raised to -_k_.
#' @param dim_names Dimnames of adjacency matrix
#'
#' @returns Column-normalized matrix.
#'
column_normalize = function(adjM, norm, k, dim_names){
  if(any(unlist(lapply(dim_names, function(x) length(x) == 1)))){
    # rL = length(dim_names[[1]])
    # cL = length(dim_names[[2]])
    # adjM = Matrix::Matrix(adjM, nrow = rL, ncol = cL, sparse = TRUE, dimnames = dim_names)
    dimnames(adjM) = dim_names
    adjM = methods::as(methods::as(methods::as(adjM, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  }
  if(norm == "degree"){
    wt <- Matrix::Diagonal(x = Matrix::colSums(adjM)^(-1))
    nadjM <- adjM %*% wt
  }else if(norm == "modified_degree"){
    adjM.mod = modify_adj_mat(adjM = adjM, k = k)
    wt.mod <- Matrix::Diagonal(x = Matrix::colSums(adjM.mod)^(-1))
    nadjM <- adjM.mod %*% wt.mod
  }
  dimnames(nadjM) = dim_names
  sum2one(nadjM)
}

#' @title Coerce columns of a matrix to sum to one
#'
#' @description
#' Divides columns by column sums so that all columns sum to one.
#'
#' @param X matrix
#'
#' @returns A matrix, with all columns summing to one.
#'
sum2one <- function(X) {
  x.dimnames = dimnames(X)
  if(!"dgCMatrix" %in% class(X)) X = methods::as(methods::as(methods::as(X, "dMatrix"), "generalMatrix"), "CsparseMatrix") # X = Matrix::Matrix(X, sparse = TRUE)
  inv_col_sum = Matrix::Diagonal(x = Matrix::colSums(X)^(-1))
  res = X %*% inv_col_sum
  dimnames(res) = x.dimnames
  res[is.na(res)] = 0 # This ensures that columns with all zeros that were divided by zero remain zero.
  res
}

#' @title Extract the diagonal elements of a matrix
#'
#' @description Extract the diagonal elements of a square matrix.
#'
#' @param X matrix
#'
#' @returns vector of diagonal elements of X.
#'
#' @export
get_diagonal = function(X){
  if(!"dgCMatrix" %in% class(X)) X = methods::as(methods::as(methods::as(X, "dMatrix"), "generalMatrix"), "CsparseMatrix")
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

#' @title Iterative Proportional Fitting
#'
#' @description Scales input matrix X to have row and column sums approximately equal to 1, thereby transforming X to a bistochastic (doubly stochastic) matrix.
#'
#' @param X Matrix
#' @param e Edge weight to add to zero-value diagonal elements of X to aid on convergence
#'
#' @returns A named list
#'  B: scaled matrix from the matrix product PXQ
#'  p: diagonal elements of diagonal matrix P
#'  q: diagonal elements of diagonal matrix Q
#'
#'  @export
ipf = function(X, e = 1e-6){
  gamma = 1
  # For starting transition matrix X, obtain B = PXQ, where B is bistochastic, P,Q are diagonal matrices, and B and X are as similar as possible (minimize relative entropy between B & X)
  # Set target row sums depending on how much to homogenize degree influence. gamma = 1 for complete homogenization, gamma = 0 for no correction.
  get_rowsum_targets = function(rs, l, t) rs - l * (rs - t)
  target.row.sums = get_rowsum_targets(Matrix::rowSums(X), gamma, 1)
  target.col.sums = 1

  x.dimnames = dimnames(X)

  # Adding self-loops to aid in convergence
  d = get_diagonal(X)
  d[d == 0] = e
  X = X + Matrix::Diagonal(n = nrow(X), x = d)

  stop_delta = 1e-6
  step = 1
  stop_step = 200
  q1 = rep(1, nrow(X)) # initialize the diagonal elements of Q
  p1 = rep(1, nrow(X)) # initialize the diagonal elements of P
  repeat{
    # set p, given q (p are the diagonal elements of P)
    p2 = target.row.sums / Matrix::rowSums(X %*% Matrix::Diagonal(x = q1))
    # set q, given p found above (q are the diagonal elements of Q)
    q2 = target.col.sums / Matrix::colSums(Matrix::Diagonal(x = p2) %*% X)

    # delta.p = all(abs(p2 - p1) <= stop_delta)
    # delta.q = all(abs(q2 - q1) <= stop_delta)
    delta.p = sum(abs(p2 - p1)) <= stop_delta
    delta.q = sum(abs(q2 - q1)) <= stop_delta
    step = step + 1
    if((delta.p && delta.q) || step > stop_step) break
    q1 = q2
    p1 = p2
  }
  P = Matrix::Diagonal(x = p2)
  Q = Matrix::Diagonal(x = q2)
  B = P %*% X %*% Q
  dimnames(B) = x.dimnames
  return(list(B = B, p = p2, q = q2))
}

#' @title Bistochastic Scaling
#'
#' @description Directly modify the transition matrix to attenuate the influence of degree on diffusion scores (as evidenced by an increased entropy of stationary distribution associated with the modified transition matrix).
#' This is done by scaling the transition matrix to be approximately bistochastic (all row & column sums equal 1).
#'
#' @param trans_mat transition matrix
#'
#' @returns A modified transition matrix that is approximately bistochastic and thus has been adjusted for degree bias.
#'
#' @export
bistochastic_scaling = function(trans_mat){
  gamma = 1
  if(gamma == 0) return(trans_mat)
  x = 1e-6
  if(!"dgCMatrix" %in% class(trans_mat)) trans_mat = methods::as(methods::as(methods::as(trans_mat, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  B.tmp = ipf(trans_mat, x) # ipf(trans_mat, x, gamma)
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
  B = Matrix::drop0(B)
  return(B)
}

#' @title Inflation-Normalization Procedure
#'
#' @description Directly modify the transition matrix to attenuate the influence of degree on diffusion scores (as evidenced by an increased entropy of stationary distribution associated with the modified transition matrix).
#' This is done by raising the values in each row of a left-stochastic transition matrix by an exponent greater than one which is a positive linear function of the stationary probability of the node associated with that row. This is followed by column normalization. This procedure displaces incoming transition probabilities to a node to the other outgoing transition probabilities of its neighbors as a function of degree.
#'
#' @param nadjM a column-normalized adjacency matrix (i.e., transition matrix)
#'
#' @returns A modified transition matrix that has been adjusted for degree bias.
#'
#' @export
inflate_normalize <- function (nadjM){
  # Changing to a Row-compressed sparse matrix
  nadjM = methods::as(methods::as(methods::as(nadjM, "dMatrix"), "generalMatrix"), "RsparseMatrix")
  # Getting stationary distribution of nadjM
  stat.distr1 = stationary.distr(nadjM)
  e0 = entropy(stat.distr1)
  # Performing Inflation/normalization on transition matrix
  kf = c(1, 10, seq(50, 2000, 50))
  res = numeric(length(kf))
  for(j in seq_along(kf)){
    inflation = 1 + kf[j] * stat.distr1
    nadjM.tmp = nadjM
    tmp = numeric(Matrix::nnzero(nadjM.tmp))
    for(i in which(Matrix::rowSums(nadjM.tmp) != 0)){ # i corresponds to rows of nadjM
      # id.tmp = (nadjM.tmp@p[i] + 1):nadjM.tmp@p[i+1] # ids of nadjM.tmp@x that are non-zero and in row i
      # n.tmp = nadjM.tmp@p[i+1] - nadjM.tmp@p[i] # number of non-zero elements in row i
      # tmp[id.tmp] = rep(inflation[i], n.tmp)
      tmp[(nadjM.tmp@p[i] + 1):nadjM.tmp@p[i+1]] = rep(inflation[i], nadjM.tmp@p[i+1] - nadjM.tmp@p[i])
    }
    nadjM.tmp@x = nadjM.tmp@x ^ tmp
    nadjM.tmp = sum2one(nadjM.tmp)
    sd.tmp = stationary.distr(nadjM.tmp)
    if(any(sd.tmp < 0 | sd.tmp > 1)){
      if(j == 1){
        nadjM = methods::as(methods::as(methods::as(nadjM, "dMatrix"), "generalMatrix"), "CsparseMatrix")
        return(nadjM)
      }else{
        j = j-1
        break
      }
    }
    res[j] = entropy(sd.tmp)
    if(j == 1){
      if(res[j] < e0){
        nadjM = methods::as(methods::as(methods::as(nadjM, "dMatrix"), "generalMatrix"), "CsparseMatrix")
        return(nadjM)
      }
    }else if(res[j] <= res[j-1]){
      break
    }
  }
  res = res[1:j]
  j = which.max(res)
  inflation = 1 + kf[j] * stat.distr1
  tmp = numeric(Matrix::nnzero(nadjM))
  for(i in which(Matrix::rowSums(nadjM) != 0)){ # i corresponds to rows of nadjM
    tmp[(nadjM@p[i] + 1):nadjM@p[i+1]] = rep(inflation[i], nadjM@p[i+1] - nadjM@p[i])
  }
  nadjM@x = nadjM@x ^ tmp
  nadjM = sum2one(nadjM)
  return(nadjM)
}
