#=================================#
# Multiplex-heterogeneous Networks
#=================================#
# Pick a few representative scenarios
# Save all objects associated with a specific scenario in a separate '.rda' file

# Scenario 1: mono-hetero as edge lists, with 'data' as list
# Scenario 2: multi-hetero as list of graphs, with biased random walk attribute and data attribute
# Scenario 3: multi-homo as list of adjacency matrices

rm(list = ls())

library(igraph)
library(data.table)

# Helper functions
extract_string = function(x, k, pos) unlist(lapply(strsplit(x, k), function(y) y[pos]))
largest_connected_component = function(g){
  if(!igraph::is_connected(g)){
    comps = igraph::components(g)
    largest_comp_id = which.max(comps$csize)
    g = igraph::induced_subgraph(g, which(comps$membership == largest_comp_id))
  }
  g
}

scenario = 1

# Scenario 1 ----
if(scenario == 1){
  # Monoplex-heterogeneous network
  # List of edge lists

  file_paths = c("/Users/samboyd/Documents/GRA/Network Analysis/AMEND/R_package/10090.protein.links.v11.0.txt.gz",
                 "/Users/samboyd/Documents/GRA/Network Analysis/AMEND/R_package/chemical_chemical.links.detailed.v5.0.tsv.gz",
                 "/Users/samboyd/Documents/GRA/Network Analysis/AMEND/R_package/10090.protein_chemical.links.v5.0.tsv.gz")
  ppi = data.table::fread(file = file_paths[1],
                          header = TRUE) %>%
    dplyr::filter(combined_score >= 800) %>%
    dplyr::rename(node1 = protein1,
                  node2 = protein2) %>%
    dplyr::mutate(node1 = paste0(node1, "|prot"),
                  node2 = paste0(node2, "|prot"),
                  combined_score = combined_score / 1000) %>%
    as.matrix()
  g.prot = igraph::graph_from_edgelist(ppi[,1:2], directed=FALSE)
  E(g.prot)$weight = as.numeric(ppi[,3])
  g.prot = igraph::simplify(graph = g.prot)
  uniq.prot = unique(extract_string(c(ppi[,1], ppi[,2]), "\\|", 1))

  mmi = data.table::fread(file = file_paths[2],
                          header = TRUE) %>%
    dplyr::filter(combined_score >= 800) %>%
    dplyr::select(chemical1, chemical2, combined_score) %>%
    dplyr::rename(node1 = chemical1,
                  node2 = chemical2) %>%
    dplyr::mutate(node1 = paste0(node1, "|meta"),
                  node2 = paste0(node2, "|meta"),
                  combined_score = combined_score / 1000) %>%
    as.matrix()
  g.meta = igraph::graph_from_edgelist(mmi[,1:2], directed=FALSE)
  E(g.meta)$weight = as.numeric(mmi[,3])
  g.meta = igraph::simplify(graph = g.meta, edge.attr.comb = list(weight = "min"))
  set.seed(55)
  clust = igraph::cluster_louvain(graph = g.meta)
  c.tbl = sort(table(clust$membership), decreasing = TRUE)
  clust.id = as.numeric(names(c.tbl)[3])
  subg.meta = igraph::induced_subgraph(graph = g.meta, vids = which(clust$membership %in% clust.id))
  if(!is_connected(subg.meta)) subg.meta = largest_connected_component(subg.meta)
  mmi = as_edgelist(subg.meta)
  mmi = cbind(mmi, E(subg.meta)$weight)
  uniq.meta = unique(extract_string(c(mmi[,1], mmi[,2]), "\\|", 1))

  bp = data.table::fread(file = file_paths[3],
                         header = TRUE) %>%
    dplyr::filter(combined_score >= 800) %>%
    dplyr::filter(protein %in% uniq.prot,
                  chemical %in% uniq.meta) %>%
    dplyr::rename(node1 = chemical,
                  node2 = protein) %>%
    dplyr::mutate(node1 = paste0(node1, "|meta"),
                  node2 = paste0(node2, "|prot"),
                  combined_score = combined_score / 1000) %>%
    as.matrix()

  el = rbind(ppi, mmi, bp)

  g = igraph::graph_from_edgelist(el = el[,1:2], directed = FALSE)
  E(g)$weight = as.numeric(el[,3])
  g = igraph::simplify(graph = g, edge.attr.comb = list(weight = "min"))
  if(!is_connected(g)) g = largest_connected_component(g)
  # table(extract_string(V(g)$name, "\\|", 2))

  # To get a smaller subgraph
  if(0){
    set.seed(55)
    clust = igraph::cluster_louvain(graph = g)
    clust.id = as.numeric(names(sort(table(clust$membership), decreasing = TRUE))[1])
    subg = igraph::induced_subgraph(graph = g, vids = which(clust$membership %in% clust.id))
    c(vcount(subg), ecount(subg))
    if(!is_connected(subg)) subg = largest_connected_component(subg)
    g = subg
    # table(extract_string(V(subg)$name, "\\|", 2))
  }

  bp.node = logical(vcount(g))
  for(i in seq_along(bp.node)){
    nei = as_ids(neighbors(g, i, mode="all"))
    bp.node[i] = all(c("prot", "meta") %in% extract_string(nei, "\\|", 2))
  }

  ppi = igraph::as_edgelist(induced_subgraph(g, grep(pattern = "\\|prot", x = V(g)$name)))
  ppi = cbind(ppi, E(induced_subgraph(g, grep(pattern = "\\|prot", x = V(g)$name)))$weight)
  mmi = igraph::as_edgelist(induced_subgraph(g, grep(pattern = "\\|meta", x = V(g)$name)))
  mmi = cbind(mmi, E(induced_subgraph(g, grep(pattern = "\\|meta", x = V(g)$name)))$weight)
  bp = igraph::as_edgelist(induced_subgraph(g, which(bp.node)))
  bp = cbind(bp, E(induced_subgraph(g, which(bp.node)))$weight)

  edgelists = list(prot = ppi, meta = mmi, "prot;meta" = bp)

  save(edgelists, file = "/Users/samboyd/Documents/GRA/Network Analysis/AMEND/R_package/AMEND/data/mono_hetero_edgelists.rda")

  # TEST
  if(0){
    rm(list=ls())
    library(igraph)
    # Create object for 'data' argument
    data.list.input = vector("list", 2)
    names(data.list.input) = c("prot", "meta")
    for(i in seq_along(data.list.input)){
      uniq.names = extract_string(unique(c(edgelists[[i]][,1], edgelists[[i]][,2])), "\\|", 1)
      # data.list.input[[i]] = rexp(n = length(uniq.names), rate = 2) # NULL
      data.list.input[[i]] = rnorm(n = length(uniq.names), mean = 0, sd = 1) # shift_scale or exp
      # data.list.input[[i]] = sample(x = c(0,1), size = length(uniq.names), replace = T, prob = c(0.9, 0.1)) # binary
      # data.list.input[[i]] = runif(length(uniq.names)) # p-value
      names(data.list.input[[i]]) = uniq.names
    }

    # Setting multiplex/heterogeneous parameters
    net.weight = c(prot = 0.5, meta = 0.5)
    jump.prob = c(prot = 0.7, meta = 0.3)

    # FUN and FUN.params
    FUN = list(prot = "shift_scale", meta = function(x, k) abs(x) * k)
    FUN.params = list(prot = list(DOI = 1, w = 0.5), meta = list(k = 1.2))

    load("/Users/samboyd/Documents/GRA/Network Analysis/AMEND/R_package/AMEND/data/mono_hetero_edgelists.rda")
    devtools::load_all()
    g = create_integrated_graph(edge_list = edgelists, data = data.list.input, node_type = NULL, brw.attr = NULL,
                                FUN = FUN, FUN.params = FUN.params, heterogeneous = TRUE, multiplex = FALSE)
    subnet = run_AMEND(edge_list = edgelists, n = 50, data = data.list.input, node_type = NULL, brw.attr = NULL,
                       FUN = FUN, FUN.params = FUN.params, heterogeneous = TRUE, multiplex = FALSE,
                       normalize = "degree", jump.prob = jump.prob, net.weight = net.weight, switch.layer.prob = switch.layer.prob, layer.weight = layer.weight, verbose = TRUE)
  }
}
# Scenario 2 ----
if(scenario == 2){
  file_paths = c("/Users/samboyd/Documents/GRA/Network Analysis/AMEND/R_package/10090.protein.links.v11.0.txt.gz",
                 "/Users/samboyd/Documents/GRA/Network Analysis/AMEND/R_package/chemical_chemical.links.detailed.v5.0.tsv.gz",
                 "/Users/samboyd/Documents/GRA/Network Analysis/AMEND/R_package/10090.protein_chemical.links.v5.0.tsv.gz")
  ppi = data.table::fread(file = file_paths[1],
                          header = TRUE) %>%
    dplyr::filter(combined_score >= 800) %>%
    dplyr::rename(node1 = protein1,
                  node2 = protein2) %>%
    dplyr::mutate(node1 = paste0(node1, "|prot"),
                  node2 = paste0(node2, "|prot"),
                  combined_score = combined_score / 1000) %>%
    as.matrix()
  g.prot = igraph::graph_from_edgelist(ppi[,1:2], directed=FALSE)
  E(g.prot)$weight = as.numeric(ppi[,3])
  g.prot = igraph::simplify(graph = g.prot)
  uniq.prot = unique(extract_string(c(ppi[,1], ppi[,2]), "\\|", 1))

  mmi = data.table::fread(file = file_paths[2],
                          header = TRUE) %>%
    dplyr::filter(combined_score >= 800) %>%
    dplyr::select(chemical1, chemical2, combined_score) %>%
    dplyr::rename(node1 = chemical1,
                  node2 = chemical2) %>%
    dplyr::mutate(node1 = paste0(node1, "|meta"),
                  node2 = paste0(node2, "|meta"),
                  combined_score = combined_score / 1000) %>%
    as.matrix()
  g.meta = igraph::graph_from_edgelist(mmi[,1:2], directed=FALSE)
  E(g.meta)$weight = as.numeric(mmi[,3])
  g.meta = igraph::simplify(graph = g.meta, edge.attr.comb = list(weight = "min"))
  set.seed(55)
  clust = igraph::cluster_louvain(graph = g.meta)
  c.tbl = sort(table(clust$membership), decreasing = TRUE)
  clust.id = as.numeric(names(c.tbl)[3])
  subg.meta = igraph::induced_subgraph(graph = g.meta, vids = which(clust$membership %in% clust.id))
  if(!is_connected(subg.meta)) subg.meta = largest_connected_component(subg.meta)
  mmi = as_edgelist(subg.meta)
  mmi = cbind(mmi, E(subg.meta)$weight)
  uniq.meta = unique(extract_string(c(mmi[,1], mmi[,2]), "\\|", 1))

  bp = data.table::fread(file = file_paths[3],
                         header = TRUE) %>%
    dplyr::filter(combined_score >= 800) %>%
    dplyr::filter(protein %in% uniq.prot,
                  chemical %in% uniq.meta) %>%
    dplyr::rename(node1 = chemical,
                  node2 = protein) %>%
    dplyr::mutate(node1 = paste0(node1, "|meta"),
                  node2 = paste0(node2, "|prot"),
                  combined_score = combined_score / 1000) %>%
    as.matrix()

  el = rbind(ppi, mmi, bp)

  g = igraph::graph_from_edgelist(el = el[,1:2], directed = FALSE)
  E(g)$weight = as.numeric(el[,3])
  g = igraph::simplify(graph = g, edge.attr.comb = list(weight = "min"))
  if(!is_connected(g)) g = largest_connected_component(g)
  # table(extract_string(V(g)$name, "\\|", 2))
  # To get a smaller subgraph
  if(0){
    set.seed(55)
    clust = igraph::cluster_louvain(graph = g)
    clust.id = as.numeric(names(sort(table(clust$membership), decreasing = TRUE))[1])
    subg = igraph::induced_subgraph(graph = g, vids = which(clust$membership %in% clust.id))
    c(vcount(subg), ecount(subg))
    if(!is_connected(subg)) subg = largest_connected_component(subg)
    g = subg
    # table(extract_string(V(subg)$name, "\\|", 2))
  }

  # Creating multiplex component
  g.prot = igraph::induced_subgraph(g, which(extract_string(V(g)$name, "\\|", 2) == "prot"))
  V(g.prot)$name = extract_string(V(g.prot)$name, "\\|", 1)
  if(!igraph::is_connected(g.prot)) g.prot = largest_connected_component(g.prot)
  g.tmp = vector("list", 3); names(g.tmp) = paste(rep("prot", 3), 1:3, sep = "_")
  # Split PPIN into layers for transcriptomic, proteomic, and phospho-proteomic data
  clusts = igraph::cluster_louvain(graph = g.prot) # Louvain method for getting clusters
  k = 1 # Include k-nearest neighbors
  for(i in 1:3){
    t.ids = which(clusts$membership == i)
    t.ids = unique(as.numeric(unlist(igraph::neighborhood(graph = g.prot, order = k, nodes = t.ids))))
    tmp = igraph::induced_subgraph(g.prot, t.ids)
    g.tmp[[i]] = largest_connected_component(tmp)
  }

  g.meta = igraph::induced_subgraph(g, grep(pattern = "\\|meta", x = V(g)$name))
  if(!igraph::is_connected(g.meta)) g.meta = largest_connected_component(g.meta)

  bp.node = logical(vcount(g))
  for(i in seq_along(bp.node)){
    nei = as_ids(neighbors(g, i, mode="all"))
    bp.node[i] = all(c("prot", "meta") %in% extract_string(nei, "\\|", 2))
  }

  g.bp = igraph::induced_subgraph(g, which(bp.node))

  list.of.graphs = c(g.tmp, meta = list(g.meta), "prot;meta" = list(g.bp))

  save(list.of.graphs, file = "/Users/samboyd/Documents/GRA/Network Analysis/AMEND/R_package/AMEND/data/multi_hetero_graphs.rda")

  # TEST
  if(0){
    rm(list=ls())
    library(igraph)
    load("/Users/samboyd/Documents/GRA/Network Analysis/AMEND/R_package/AMEND/data/multi_hetero_graphs.rda")
    # Create data vertex attributes, called 'scores'
    v.attr.name = "scores"
    node.types = c("prot_1", "prot_2", "prot_3", "meta")
    for(i in seq_along(node.types)){
      igraph::vertex_attr(list.of.graphs[[node.types[i]]], v.attr.name) = runif(vcount(list.of.graphs[[node.types[i]]])) # p-value
      # igraph::vertex_attr(V(list.of.graphs[[node.types[i]]]), v.attr.name) = rexp(n = length(uniq.names), rate = 2) # NULL
      # igraph::vertex_attr(V(list.of.graphs[[node.types[i]]]), v.attr.name) = rnorm(n = length(uniq.names), mean = 0, sd = 1) # shift_scale or exp
      # igraph::vertex_attr(V(list.of.graphs[[node.types[i]]]), v.attr.name) = sample(x = c(0,1), size = length(uniq.names), replace = T, prob = c(0.9, 0.1)) # binary
    }

    # Setting multiplex/heterogeneous parameters
    layer.weight = list(prot = rep(1/3, 3))
    switch.layer.prob = list(prot = c(prot_1 = 0.5, prot_2 = 0.2, prot_3 = 0.8))
    net.weight = c(prot = 0.5, meta = 0.5)
    jump.prob = c(prot = 0.5, meta = 0.3)

    devtools::load_all()
    g = create_integrated_graph(graph = list.of.graphs, data = "scores", node_type = NULL, brw.attr = NULL,
                                FUN = FUN, FUN.params = FUN.params, heterogeneous = TRUE, multiplex = TRUE)
    subnet = run_AMEND(graph = list.of.graphs, n = 100, data = "scores", node_type = NULL, brw.attr = NULL,
                       FUN = "p_value", FUN.params = NULL, heterogeneous = TRUE, multiplex = TRUE,
                       normalize = "modified_degree", k = 0.5, jump.prob = jump.prob, net.weight = net.weight,
                       switch.layer.prob = switch.layer.prob, layer.weight = layer.weight, verbose = TRUE)
  }
}
# Scenario 3 ----
if(scenario == 3){
  file_paths = c("/Users/samboyd/Documents/GRA/Network Analysis/AMEND/R_package/10090.protein.links.v11.0.txt.gz")
  ppi = data.table::fread(file = file_paths[1],
                          header = TRUE) %>%
    dplyr::filter(combined_score >= 800) %>%
    dplyr::rename(node1 = protein1,
                  node2 = protein2) %>%
    dplyr::mutate(combined_score = combined_score / 1000) %>%
    as.matrix()
  g.prot = igraph::graph_from_edgelist(ppi[,1:2], directed=FALSE)
  E(g.prot)$weight = as.numeric(ppi[,3])
  g.prot = igraph::simplify(graph = g.prot, edge.attr.comb = list(weight = "min"))
  if(!igraph::is_connected(g.prot)) g.prot = largest_connected_component(g.prot)

  # Creating multiplex component
  set.seed(305)
  list.of.adjmats = vector("list", 3); names(list.of.adjmats) = paste(rep("prot", 3), 1:3, sep = "_")
  # Split PPIN into layers for transcriptomic, proteomic, and phospho-proteomic data
  clusts = igraph::cluster_louvain(graph = g.prot) # Louvain method for getting clusters
  k = 1 # Include k-nearest neighbors
  for(i in 1:3){
    t.ids = which(clusts$membership == i)
    t.ids = unique(as.numeric(unlist(igraph::neighborhood(graph = g.prot, order = k, nodes = t.ids))))
    tmp = igraph::induced_subgraph(g.prot, t.ids)
    tmp = largest_connected_component(tmp)
    list.of.adjmats[[i]] = igraph::as_adjacency_matrix(tmp, attr = "weight", sparse = T)
  }

  save(list.of.adjmats, file = "/Users/samboyd/Documents/GRA/Network Analysis/AMEND/R_package/AMEND/data/multi_homo_adjmats.rda")

  # TEST
  if(0){
    rm(list=ls())
    library(igraph)
    load("/Users/samboyd/Documents/GRA/Network Analysis/AMEND/R_package/AMEND/data/multi_hetero_graphs.rda")

    # Create object for 'data' argument
    data.list.input = vector("list", 3)
    names(data.list.input) = c("prot_1", "prot_2", "prot_3")
    for(i in seq_along(data.list.input)){
      # data.list.input[[i]] = rexp(n = nrow(list.of.adjmats[[i]]), rate = 2) # NULL
      # data.list.input[[i]] = rnorm(n = nrow(list.of.adjmats[[i]]), mean = 0, sd = 1) # shift_scale or exp
      data.list.input[[i]] = sample(x = c(0,1), size = nrow(list.of.adjmats[[i]]), replace = T, prob = c(0.9, 0.1)) # binary
      # data.list.input[[i]] = runif(nrow(list.of.adjmats[[i]])) # p-value
      names(data.list.input[[i]]) = rownames(list.of.adjmats[[i]])
    }

    # Setting multiplex/heterogeneous parameters
    layer.weight = list(prot = c(prot_1 = 0.4, prot_2 = 0.3, prot_3 = 0.3))
    switch.layer.prob = list(prot = c(prot_1 = 0.4, prot_2 = 0.5, prot_3 = 0.5))

    devtools::load_all()
    g = create_integrated_graph(adj_matrix = list.of.adjmats, data = data.list.input, node_type = NULL, brw.attr = NULL,
                                FUN = "binary", FUN.params = NULL, heterogeneous = FALSE, multiplex = TRUE)
    subnet = run_AMEND(adj_matrix = list.of.adjmats, n = 100, data = data.list.input, node_type = NULL, brw.attr = NULL,
                       FUN = "binary", FUN.params = NULL, heterogeneous = FALSE, multiplex = TRUE, normalize = "modified_degree", k = 0.5,
                       switch.layer.prob = switch.layer.prob, layer.weight = layer.weight, verbose = TRUE)
  }
}



