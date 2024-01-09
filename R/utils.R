#' @title Extract substrings before or after a delimiter
#'
#' @description
#' This facilitates identifying the component and the layer where each node is located by splitting on '|' and '`_`' delimiters. This component/layer information is appended to node names in the form '|component_layer'.
#'
#' @param x vector of character strings
#' @param k Pattern to split on
#' @param before Return everything before the split or after the split?
#'
#' @returns character string
#'
extract_string = function(x, k, before = T) unlist(lapply(strsplit(x, k), function(y) y[ifelse(before,1,2)]))
# extract_string = function(x, k, pos) unlist(lapply(strsplit(x, k), function(y) y[pos]))

#' @title Get network location info from node labels
#'
#' @description
#' Retrieves either the node name (_mode_=1), the node component and layer (_mode_=2), just the component (_mode_=3), or just the layer (_mode_=4).
#'
#' Node labels follow the form 'name|component_layer'.
#'
#' @param x Node labels.
#' @param mode Integer. Type of node info to return. 1=node name. 2=component and layer names of nodes in _x_. 3=component names of nodes in _x_. 4=layer names of nodes in _x_.
#'
#' @returns Character vector of info specified by _mode_, extracted from node labels.
#'
get.type = function(x, mode){
  if(mode == 1){ # get node name
    t = extract_string(x, "\\|", before=TRUE)
    # t = extract_string(x, "\\|", pos = 1)
  }
  if(mode == 2){ # get node component_layer
    t = extract_string(x, "\\|", before=FALSE)
    # t = extract_string(x, "\\|", pos = 2)
  }
  if(mode == 3){ # get node component
    t = extract_string(extract_string(x, "\\|", before=FALSE), "_", before=TRUE)
    # t = extract_string(extract_string(x, "\\|", pos = 2), "_", pos = 1)
  }
  if(mode == 4){ # get node layer
    t = extract_string(extract_string(x, "\\|", before=FALSE), "_", before=FALSE)
    # t = extract_string(extract_string(x, "\\|", pos = 2), "_", pos = 2)
  }
  t
}
