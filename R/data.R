#' @title Mus musculus PPI network
#'
#' @description
#' A subset of the Mus musculus protein-protein interaction network, including both functional and physical interactions. Only interactions with a combined score greater than or equal to 0.8 were included
#'
#' @format an igraph object containing an edge attribute "weight" and the following vertex attributes:
#' * name: Ensmebl peptide ID
#' * symbol: MGI symbols
#' * ECI: Equivalent change indices (ECI)
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35378}, \url{https://version-11-0b.string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Mus+musculus}
"glut4_graph"


#' @title Mus musculus PPI network adjacency matrix
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35378}, \url{https://version-11-0b.string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Mus+musculus}
"glut4_adjM"

#' @title ECI scores for proteins in the PPI network
#'
#' @description
#' Equivalent Change Index (ECI) scores from a GLUT4 knockout-overexpression (KO-OX) microarray experiment in mice. The ECIs are ratio's of log fold changes from a differential expression analysis, in this case from the KO-Ctrl and OX-Ctrl comparisons
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35378}, \url{https://version-11-0b.string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Mus+musculus}
"eci_scores"

#' @title Log fold changes for GLUT4 KO vs. Control for proteins in the PPI network
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35378}, \url{https://version-11-0b.string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Mus+musculus}
"logFC_KO"

#' @title List of edge lists representing a monoplex-heterogeneous network
#'
#' @description
#' Protein-protein, chemical-chemical, and protein-chemical interactions, along with interaction confidence scores, in the form of edge lists. An edge list is a 2-column matrix with each row denoting an interaction. A possible third row can contain interaction scores.
#'
#' @format
#' A named list with 3 elements:
#'  * 'prot': protein-protein interaction edge list
#'  * 'meta': chemical-chemical interaction edge list (meta=metabolite=chemical)
#'  * 'prot;meta': protein-chemical interaction edge list
#'
#' @source
#' Protein-protein interactions come from STRING. Chemical-chemical and chemical-protein interactions come from STITCH.
#'
"edgelists"

#' @title List of igraphs representing a multiplex-heterogeneous network
#'
#' @description
#' 5 separate igraph objects representing layers of a multiplex protein-protein interaction network, a chemical-chemical interaction network, and a bipartite chemical-protein interaction network. The multiplex layers were taken from 3 different clusters of the full PPI network after running the Louvain clustering algorithm.
#'
#' @format
#' A named list with 5 elements:
#'  * 'prot_1': igraph of layer 1 of protein multiplex component
#'  * 'prot_2': igraph of layer 2 of protein multiplex component
#'  * 'prot_3': igraph of layer 3 of protein multiplex component
#'  * 'meta': chemical-chemical interaction igraph
#'  * 'prot;meta': igraph of bipartite protein-chemical interaction network
#'
#' @source
#' Protein-protein interactions come from STRING. Chemical-chemical and chemical-protein interactions come from STITCH.
#'
"list.of.graphs"

#' @title List of adjacency matrices representing a multiplex-homogeneous network
#'
#' @description
#' 3 adjacency matrices representing 3 layers of a multiplex protein-protein interaction network. The multiplex layers were taken from 3 different clusters of the full PPI network after running the Louvain clustering algorithm.
#'
#' @format
#' A named list with 3 elements:
#'  * 'prot_1': adjacency matrix of layer 1 of protein multiplex component
#'  * 'prot_2': adjacency matrix of layer 2 of protein multiplex component
#'  * 'prot_3': adjacency matrix of layer 3 of protein multiplex component
#'
#' @source
#' Protein-protein interactions come from STRING.
#'
"list.of.adjmats"
