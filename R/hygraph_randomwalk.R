#' hygraph_randomwalk
#' 
#' hygraph_randomwalk calculates the importance score for each sample using pre-processed data
#' 
#' @param sample_i the i-th sample
#' @param data_mut_com usable mutation data
#' @param com_mut other samples that have co-mutated gene(s) with the i-th sample
#' @param correlation_samples similarity matrix for the samples
#' @param Graph the graph object in igraph format
#' 
#' @return importance score for i-th sample
#' @importFrom utils read.delim
#' @examples
#' 
#' 
#' @export hygraph_randomwalk
hygraph_randomwalk<-
  function(sample_i,data_mut_com,co_mut,correlation_samples,Graph){
  
  sample_i_mut_genes=rownames(data_mut_com)[which(data_mut_com[,sample_i]==1)] 
  
  com_mut_sample_i=co_mut[sample_i,which(co_mut[sample_i,]!=0)]  
  com_mut_samples=names(com_mut_sample_i)
  
  H=data_mut_com[sample_i_mut_genes,com_mut_samples,drop=F] 
  #dim(H)
  Vertex_W=get_vertex_W(sample_i,data_mut_com,com_mut_samples,Graph,sample_i_mut_genes,H)
  
  hyperedge_W_vector=correlation_samples[,sample_i]
  hyperedge_W=diag(hyperedge_W_vector[colnames(H)])
  colnames(hyperedge_W)=colnames(H)
  rownames(hyperedge_W)=colnames(H)
  
  P=get_P(H,Vertex_W,hyperedge_W)
  importance_score = get_hyper_randomwalk(P,H,theta=0.2)
  
  return(importance_score)
}
