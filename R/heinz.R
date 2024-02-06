#' @importFrom igraph vcount V V<- E E<- vertex_attr vertex_attr<-

#' @title Heuristic solution of Maximum-weight Connected Subgraph problem
#'
#' @description Given a graph and a named vector of node scores, `heinz()` heuristically finds a solution to the maximum-weight connected subgraph problem.
#'
#' @details Details can be found in the _dnet_ package documentation under the `dnetfind()` function: \url{https://cran.r-project.org/web/packages/dnet/dnet.pdf}
#'
#' @note Adapted from the _dnet_ package. Based on method from the _BioNet_ package
#'
#' @param ig Input graph
#' @param scores Named vector of scores for nodes in ig
#' @param min.cluster.size Minimum size of disconnected components of linker graph to consider when searching for maximum-weight subgraph. If NULL, only considers largest connected component of linker graph.
#'
#' @return a subnetwork as an igraph object
#'
#' @seealso [run_AMEND()], [RWR()]
#'
#' @examples
#' library(igraph)
#'
#' graph = sample_pa(n = 100, power = 1.2)
#' V(graph)$name = 1:100
#' g_scores = rnorm(100)
#' names(g_scores) = 1:100
#'
#' new_g = heinz(ig = graph, scores = g_scores, min.cluster.size = 2)
#'
#' @export
heinz = function(ig, scores, min.cluster.size = 2){
  if(!igraph::is_igraph(ig)){
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
  # Avoiding possible vertex attribute naming conflicts
  van = "heinz_score"
  while(van %in% igraph::vertex_attr_names(ig)){
    van = paste0(van, 1)
  }
  igraph::vertex_attr(ig, van) = scores[V(ig)$name]
  pos.nodes = which(igraph::vertex_attr(ig, van) > 0)
  if(length(pos.nodes) == 0){
    warning("No positive nodes")
    subgraph <- igraph::graph.empty(n = 0, directed = FALSE)
  }else if(length(pos.nodes) == 1){
    subgraph = igraph::induced_subgraph(ig, pos.nodes)
    V(subgraph)$type <- "desired"
  }else{
    # Get the subgraphs consisting of only positive nodes
    pos.subgraph = igraph::induced_subgraph(ig, pos.nodes)
    # Returns list of connected components of positive subgraph above. These will be the meta-nodes
    conn.comp.graph <- igraph::decompose(pos.subgraph)
    # Returns the sum of the scores for each meta-node
    score.comp <- unlist(lapply(lapply(conn.comp.graph, function(x) as.numeric(igraph::vertex_attr(x, van))), sum))
    ind_order <- order(score.comp, decreasing = TRUE)
    conn.comp.graph <- conn.comp.graph[ind_order]
    score.comp <- score.comp[ind_order]
    for(i in seq_along(conn.comp.graph)){
      conn.comp.graph[[i]]$score <- score.comp[i]
    }
    v.id <- seq(1, vcount(ig))
    names(v.id) <- V(ig)$name
    edgelist <- igraph::get.edgelist(ig, names = FALSE)
    edgelist1 <- edgelist[, 1]
    edgelist2 <- edgelist[, 2]
    # This for loop is changing the edgelist to treat the nodes in each meta-node as one node
    ig.size <- vcount(ig)
    conn.comp.gt1.ids = which(unlist(lapply(conn.comp.graph, vcount)) > 1)
    conn.comp.eq1.ids = which(unlist(lapply(conn.comp.graph, vcount)) == 1)
    names(conn.comp.eq1.ids) = unlist(lapply(conn.comp.graph[conn.comp.eq1.ids], function(x) V(x)$name))
    for(i in seq_along(conn.comp.gt1.ids)){ # for each meta-node i
      # change the source, target ids of all edges containing nodes in meta-node i to be the same
      v.id.tmp = v.id[V(conn.comp.graph[[conn.comp.gt1.ids[i]]])$name]
      edgelist1[edgelist1 %in% v.id.tmp] = ig.size + i # new id
      edgelist2[edgelist2 %in% v.id.tmp] = ig.size + i # new id
    }
    # Order of new.ids and new.names need to correspond to order of conn.comp.graph
    # IDs of nodes in graph that are meta nodes connected only to negative nodes.
    v.meta.id = v.id[match(names(conn.comp.eq1.ids), names(v.id))] # v.meta.id is in same order as conn.comp.eq1.ids
    new.ids = numeric(length(conn.comp.graph))
    new.ids[c(conn.comp.eq1.ids, conn.comp.gt1.ids)] <- c(unname(v.meta.id), seq(vcount(ig) + 1, vcount(ig) + length(conn.comp.gt1.ids)))
    # For the names cluster'x', the 'x' must correspond to the xth element of conn.comp.graph
    new.names = paste("cluster", 1:length(conn.comp.graph), sep = "")
    names(new.ids) <- new.names
    v.id <- c(v.id[!v.id %in% v.meta.id], new.ids)
    v.name <- names(v.id)
    names(v.name) <- v.id
    # This is extracting from v.name only those nodes whose IDs appear in edgelist.
    new.edgelist <- cbind(v.name[as.character(edgelist1)],
                          v.name[as.character(edgelist2)])
    mig <- igraph::graph_from_edgelist(new.edgelist, directed = FALSE)
    mig <- igraph::simplify(mig, remove.loops = TRUE, remove.multiple = TRUE)
    node.score <- scores[V(mig)$name]
    names(node.score) <- V(mig)$name
    node.score.cluster <- sapply(conn.comp.graph, igraph::get.graph.attribute, "score")
    names(node.score.cluster) <- new.names
    ind_cluster <- grep("cluster", names(node.score))
    node.score[ind_cluster] <- node.score.cluster[names(node.score[ind_cluster])]
    # Appending edge weight attribute: absolute value of sum of degree-normalized node scores (only including negative nodes in this sum). Larger score = less desirable
    igraph::vertex_attr(mig, van) <- node.score
    score.degree <- 1/(igraph::degree(mig) + 1)
    tmp_score <- igraph::vertex_attr(mig, van)
    tmp_score[tmp_score > 0] <- 0
    V(mig)$score.degree <- score.degree * tmp_score
    E(mig)$weight <- rep(0, length(E(mig)))
    tmp_n1 <- igraph::get.edgelist(mig, names = FALSE)[, 1]
    tmp_n2 <- igraph::get.edgelist(mig, names = FALSE)[, 2]
    E(mig)$weight <- -(V(mig)[tmp_n1]$score.degree + V(mig)[tmp_n2]$score.degree)
    # This is to ensure the transformed graph is connected. If input graph is connected, this will always be connected.
    # If not connected, choose component with largest sum of meta-node scores
    if(!igraph::is_connected(mig)){
      decomp.graphs <- igraph::decompose(mig)
      sum.pos <- lapply(decomp.graphs, function(x) {
        sum(node.score[names(which(node.score[V(x)$name] > 0))])
      })
      mig <- decomp.graphs[[which.max(sum.pos)]]
    }
    # Minimum spanning tree of transformed graph
    mst <- igraph::minimum.spanning.tree(mig, weights = E(mig)$weight)
    mst.cluster.id <- grep("cluster", V(mst)$name)
    names(mst.cluster.id) <- V(mst)[mst.cluster.id]$name
    tmp <- unlist(strsplit(names(mst.cluster.id), "cluster"))
    ttmp <- as.numeric(matrix(tmp, nrow = 2)[2, ])
    mst.cluster.id <- mst.cluster.id[order(ttmp)]
    if(length(mst.cluster.id) == 1){
      neg.node.ids.2 = c()
    }else{
      # Iteratively exclude all negative nodes of degree 1 from MST
      loop = TRUE
      while(loop){
        neg_deg1 = which(igraph::degree(mst) == 1 & igraph::vertex_attr(mst, van) < 0)
        if(length(neg_deg1) == 0){
          loop = FALSE
        }else{
          mst = igraph::induced_subgraph(mst, -neg_deg1)
        }
      }
      # Get induced subgraph from transformed graph with these negative leaf nodes removed
      sub.mig = igraph::induced_subgraph(mig, which(V(mig)$name %in% V(mst)$name))
      if(!igraph::is_simple(sub.mig)) sub.mig = igraph::simplify(sub.mig)
      # Identify linker nodes: negative nodes that have meta-node neighbors and an absolute score less than sum of meta-node neighbor scores (Step v)
      neg.node.ids <- which(igraph::vertex_attr(sub.mig, van) < 0) # negative nodes in induced subgraph of step iv
      for(i in neg.node.ids){
        tmp_nei <- igraph::neighbors(sub.mig, v = i)
        tmp_nei_meta <- grep("cluster", V(sub.mig)[tmp_nei]$name) # If "cluster" is in the name, this means it's a meta node
        V(sub.mig)[i]$clusters <- list(tmp_nei[tmp_nei_meta])
      }
      score.neg.nodes <- c()
      for(i in neg.node.ids){
        if(!is.na(V(sub.mig)[i]$clusters[1])){ # If neg.node i has at least one meta-node neighbor...
          borders <- c(i, V(sub.mig)[i]$clusters) # vector of id of neg.node and ids of its meta-node neighbors
          borders <- unlist(borders)
          score.neg.nodes = c(score.neg.nodes, sum(igraph::vertex_attr(sub.mig, van, borders))) # sum of scores of neg.node and its meta-node neighbors
        }else{
          score.neg.nodes <- c(score.neg.nodes, igraph::vertex_attr(sub.mig, van, i))
        }
      }
      neg.node.ids.2 <- neg.node.ids[score.neg.nodes >= 0] # Or neg.node.ids[score.neg.nodes >= 0], to favor larger subnetworks in case of ties
      # If score.neg.nodes is negative for a neg-node, this necessarily means it wasn't connected to any meta nodes and will be removed.
    }
    # Getting MST of linker graph (containing linker nodes only) (step vi)
    if(length(neg.node.ids.2) == 0){
      tmp <- unlist(strsplit(names(node.score.cluster)[which.max(node.score.cluster)], "cluster"))
      ttmp <- as.numeric(matrix(tmp, nrow = 2)[2, ])
      tmp_nodes <- unlist(lapply(conn.comp.graph, igraph::get.vertex.attribute, "name")[ttmp])
      subgraph = igraph::induced_subgraph(ig, which(V(ig)$name %in% tmp_nodes))
      subgraph = largest_connected_component(subgraph)
      V(subgraph)$type <- "desired"
      subgraph = igraph::delete_vertex_attr(subgraph, van)
      return(subgraph)
    }else{
      subg = igraph::induced_subgraph(sub.mig, neg.node.ids.2) # 'linker' graph
      if(!igraph::is_connected(subg)){
        clust <- igraph::components(subg)
        max.cluster.size = max(clust$csize)
        if(min.cluster.size > max.cluster.size || is.null(min.cluster.size)) min.cluster.size = max.cluster.size
        clust.ids = which(clust$csize >= min.cluster.size)
        getPathScore <- function(path, g1.score, g1.clust, g2.score){
          s1 <- g1.score[path] # scores from a path in the MST linker graph
          tmp <- unique(unlist(g1.clust[path])) # meta nodes that linker nodes along 'path' are connected to
          s2 <- g2.score[tmp] # scores of meta nodes that linker nodes along 'path' are connected to
          sum(s1, s2)
        }
        subgraph.list = vector("list", length(clust.ids))
        for(j in seq_along(clust.ids)){
          subg.tmp <- igraph::induced_subgraph(subg, which(clust$membership == clust.ids[j]))
          if(!igraph::is_simple(subg.tmp)) subg.tmp = igraph::simplify(subg.tmp)
          mst.subg.tmp <- igraph::minimum.spanning.tree(subg.tmp, E(subg.tmp)$weight) # MST of linker graph
          # Find highest scoring path (by vertex weights) in the linker MST plus attached meta-nodes (Step vii in dNetFind function documentation)
          max.score <- 0
          mst.subg.score.tmp = igraph::vertex_attr(mst.subg.tmp, van)
          mst.subg.clust.tmp = V(mst.subg.tmp)$clusters
          sub.mig.score.tmp = igraph::vertex_attr(sub.mig, van)
          for(i in 1:vcount(mst.subg.tmp)){
            path <- igraph::all_shortest_paths(mst.subg.tmp, from = V(mst.subg.tmp)[i])
            path.score <- unlist(lapply(path$res, getPathScore, g1.score = mst.subg.score.tmp, g1.clust = mst.subg.clust.tmp, g2.score = sub.mig.score.tmp))
            best.pos <- which.max(path.score)
            if(path.score[[best.pos]] > max.score){
              best.path <- path$res[[best.pos]]
              max.score <- path.score[[best.pos]]
            }
          }
          if(length(best.path) != 1){
            cluster.list <- V(mst.subg.tmp)[best.path]$clusters
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
          # Induce a subgraph from input graph containing only
          # nodes along this best path and edges between them
          pos.cluster <- V(sub.mig)[unique(unlist(V(mst.subg.tmp)[best.path]$clusters))]$name
          tmp <- unlist(strsplit(pos.cluster, "cluster"))
          ttmp <- as.numeric(matrix(tmp, nrow = 2)[2, ])
          tmp_meta_nodes <- unlist(lapply(conn.comp.graph, igraph::get.vertex.attribute, "name")[ttmp])
          tmp_border_nodes <- V(mst.subg.tmp)[best.path]$name
          tmp_nodes <- c(tmp_border_nodes, tmp_meta_nodes)
          subgraph <- igraph::induced_subgraph(ig, vids = tmp_nodes)
          if(!igraph::is_connected(subgraph)){
            clust.tmp <- igraph::components(subgraph)
            cid.tmp <- which.max(clust.tmp$csize)
            subgraph <- igraph::induced_subgraph(subgraph, which(clust.tmp$membership == cid.tmp))
          }
          type <- rep("desired", vcount(subgraph))
          names(type) <- V(subgraph)$name
          type[tmp_border_nodes[tmp_border_nodes %in% names(type)]] <- "linker"
          V(subgraph)$type <- type
          subgraph.list[[j]] = subgraph
        }
        subgraph.list = subgraph.list[order(unlist(lapply(subgraph.list, vcount)), decreasing = TRUE)] # graphs are in descending order
        subgraph.scores = lapply(subgraph.list, function(x) sum(igraph::vertex_attr(x, van)))
        g.tmp = subgraph.list[[1]]
        score.tmp = sum(igraph::vertex_attr(subgraph.list[[1]], van))
        ll = ifelse(length(subgraph.list) > 1, 2, 1)
        for(l in ll:length(subgraph.list)){
          all.nodes = unique(c(V(g.tmp)$name, V(subgraph.list[[l]])$name))
          g.tmp2 = igraph::induced_subgraph(ig, which(V(ig)$name %in% all.nodes))
          score.tmp2 = sum(igraph::vertex_attr(g.tmp2, van))
          if(score.tmp2 >= score.tmp){
            g.tmp = g.tmp2
            score.tmp = score.tmp2
          }
        }
        g.tmp = igraph::delete_vertex_attr(g.tmp, van)
        g.tmp = largest_connected_component(g.tmp)
        return(g.tmp)
      }else{
        if(!igraph::is_simple(subg)) subg = igraph::simplify(subg)
        mst.subg <- igraph::minimum.spanning.tree(subg, E(subg)$weight) # MST of linker graph
        getPathScore <- function(path, g1.score, g1.clust, g2.score){
          s1 <- g1.score[path] # scores from a path in the MST linker graph
          tmp <- unique(unlist(g1.clust[path])) # meta nodes that linker nodes along 'path' are connected to
          s2 <- g2.score[tmp] # scores of meta nodes that linker nodes along 'path' are connected to
          sum(s1, s2)
        }
        # Find highest scoring path (by vertex weights) in the linker MST plus attached meta-nodes (Step vii in dNetFind function documentation)
        max.score <- 0
        mst.subg.score.tmp = igraph::vertex_attr(mst.subg, van)
        mst.subg.clust.tmp = V(mst.subg)$clusters
        sub.mig.score.tmp = igraph::vertex_attr(sub.mig, van)
        for(i in 1:vcount(mst.subg)){
          path <- igraph::all_shortest_paths(mst.subg, from = V(mst.subg)[i])
          path.score <- unlist(lapply(path$res, getPathScore, g1.score = mst.subg.score.tmp, g1.clust = mst.subg.clust.tmp, g2.score = sub.mig.score.tmp))
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
        # Induce a subgraph from input graph containing only
        # nodes along this best path and edges between them
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
        subgraph = igraph::delete_vertex_attr(subgraph, van)
        subgraph = largest_connected_component(subgraph)
        return(subgraph)
      }
    }
  }
}
