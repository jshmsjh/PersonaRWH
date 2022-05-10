#' get_vertex_W
#' 
#' get_vertex_W uses the mutation data and the graph constructed from network data to calculate the weight of the hypergraph
#' output of get_vertex_W is the vertex weight of the hypergraph
#' 
#' @param sample_i the i-th sample
#' @param data_mut_com usable mutation data
#' @param com_mut_samples other samples that have co-mutated gene(s) with the i-th sample
#' @param Graph the graph object in igraph format
#' @param sample_i_mut_genes placeholder
#' @param H adjacent matrix for i-th sample
#' @importFrom igraph V
#' @importFrom igraph induced_subgraph
#' @importFrom igraph degree
#' @importFrom igraph graph.data.frame
#' @examples
#' 
#' 
#' @export get_vertex_W
get_vertex_W <-
function(sample_i,data_mut_com,com_mut_samples,Graph,sample_i_mut_genes,H){
  
  mut_data_i=data_mut_com[,sample_i] 
  vertex_W=c()
  for(k in com_mut_samples){
    
    vertex_w=rep(0,length(sample_i_mut_genes))
    names(vertex_w)=sample_i_mut_genes
    
    
    mut_data_k=data_mut_com[,k]
    mut_gene=names(mut_data_k[which(mut_data_k==1)])
    
    
    v_i=intersect(mut_gene,V(Graph)$name)
    v_diff=setdiff(mut_gene,V(Graph)$name) 
    
    graph_i=induced_subgraph(Graph, v_i) 
    graph_i_degree=degree(graph_i)[v_i] 
    graph_i_degree[which(graph_i_degree==0)]=0.01 
    
    
    co_mut_genes=intersect(v_i,sample_i_mut_genes)
    
    
    vertex_w[co_mut_genes]=graph_i_degree[co_mut_genes]
    
    vertex_w[setdiff(sample_i_mut_genes,co_mut_genes)]=0.01
    vertex_W=cbind(vertex_W,vertex_w)
  }
  
  colnames(vertex_W)=com_mut_samples
  return(vertex_W)
}
