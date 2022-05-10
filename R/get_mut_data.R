#' get_mut_data
#' 
#' get_mut_data uses reads the mutation file
#' it deletes some row(s) of the mutation data whose rowsum is zero, rownames of the remaining rows are considered possibale genes
#' 
#' @param mutation_file the text that contains the mutation data, a gene*sample matrix with rownames being gene names and column names being sample IDs
#' 
#' @return Mutation, the mutation data, rows represent genes and columns represent samples
#' @importFrom utils read.delim
#' @examples

#' 
#' 
#' @export get_mut_data
get_mut_data <-
function(mutation_file){

  data=read.delim(mutation_file,header = T,as.is = T,check.names = F)

  especial_gene_index=which(is.na(data[,1])|data[,1]=="")
  if (length(especial_gene_index)>0){
    data=data[-especial_gene_index,]
  } else {data=data }

  rownames(data)=data[,1]
  data=data[,-1]

  mut_data=as.matrix(data[which(apply(data,1,sum)!=0),which(apply(data,2,sum)!=0)])
  
  return(mut_data)
}
