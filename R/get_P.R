#' get_P
#' 
#' get_P uses the mutation data, expression data and the vertex weight of the hypergraph to calculate the transition possibility matrix
#' output of get_P is the transition possibility matrix of the hypergraph
#' 
#' @param H adjacent matrix for i-th sample
#' @param Mutation the mutation data, rows represent genes and columns represent samples
#' @param vertex_W the vertex weight of the hypergraph
#' @examples
#' 
#' @export get_P
get_P <-
  function(H,vertex_W,hyperedge_W){
    
    degree_H = H%*%hyperedge_W
    
    Degree_v=apply(degree_H,1,sum)
    D_v_inverse=diag(1/Degree_v)
    
    probability_hyperedge=D_v_inverse%*%H%*%hyperedge_W
    rownames(probability_hyperedge)=rownames(degree_H)
    
    
    Degree_ve=apply(vertex_W,2,sum)
    D_ve_inverse=diag(1/Degree_ve)
    
    probability_vertex=D_ve_inverse%*%t(vertex_W)
    rownames(probability_vertex)=colnames(H)
    
    P = probability_hyperedge%*%probability_vertex
    return(as.matrix(P))
  }
