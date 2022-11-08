#' Mus musculus PPI network
#'
#' A subset of the Mus musculus protein-protein interaction network, including both functional and physical interactions. Only interactions with a combined score greater than or equal to 0.8 were included
#'
#' @format an igraph object containing an edge attribute "weight" and the following vertex attributes:
#' * name: Ensmebl peptide ID
#' * symbol: MGI symbols
#' * ECI: Equivalent change indices (ECI)
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35378}, \url{https://version-11-0b.string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Mus+musculus}
"glut4_graph"


#' Mus musculus PPI network adjacency matrix
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35378}, \url{https://version-11-0b.string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Mus+musculus}
"glut4_adjM"

#' ECI scores for proteins in the PPI network
#'
#' Equivalent Change Index (ECI) scores from a GLUT4 knockout-overexpression (KO-OX) microarray experiment in mice. The ECIs are ratio's of log fold changes from a differential expression analysis, in this case from the KO-Ctrl and OX-Ctrl comparisons
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35378}, \url{https://version-11-0b.string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Mus+musculus}
"eci_scores"
