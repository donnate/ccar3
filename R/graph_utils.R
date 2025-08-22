

#' Return the edge incidence matrix of an igraph graph
#'
#'
#' @param g igraph graph object.
#' @param weight edge weights.
#' 
#' @export
get_edge_incidence <- function(g, weight = 1){
  n_nodes = igraph::vcount(g)
  d_max = max(igraph::degree(g))
  #d_max = 1
  edges = data.frame(igraph::as_edgelist(g)) %>%
    dplyr::arrange(X1, X2)
  Gamma = matrix(0, nrow(edges), n_nodes)
  
  # Make beta_v into a matrix
  names_st = unique(c(edges$X1, edges$X2))
  for (e in 1:nrow(edges)){
    ind1 = which( edges$X1[e] == names_st)
    ind2 = which( edges$X2[e] == names_st)
    Gamma[e, ind1] = weight
    Gamma[e, ind2] = - weight
  }
  return(Gamma)
}