\name{PersonaRWH-package}
\alias{PersonaRWH-package}
\alias{PersonaRWH}
\docType{package}
\title{
PersonaRWH
}
\description{
PersonaRWH prioritizes altered genes by effective integration of somatic mutation data, gene expression RNAseq and molecular interaction data. A weighted Hypergraph model is formed by combining the degree of the PPI sub-network composed of each patient's mutated genes,
and the weight of hyperedges is defined by the similarity between samples. Mutated genes are ranked by the random walk on the weighted Hypergraph model for each sample.
}
\details{
\tabular{ll}{
Package: \tab PersonaRWH\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2022-05-10\cr
License: \tab MIT\cr
Authors: \tab Naiqian Zhang, Yuxuan Pang, Junhan Shi, Chenye Wang \cr
Maintainer: \tab Junhan Shi \cr
}
}
\keyword{ package }
\examples{
### the example mutation and expression data is a Lung Squamous Cell Carcinoma Dataset from TCGA
### the example network data is HumanNet
### We will show how to get example DriverRWH Results
library(PersonaRWH)

# load mutation and expression data
# you can also use get_mut_data(file) function to process mutation file in txt format
data(luscExampleMutation)
data(luscExampleExpression)

mut_data=luscExampleMutation
expr_data=luscExampleExpression

# keep the samples that have more than 2 mutation genes
mut_data_morethan_2=colnames(mut_data)[which(apply(mut_data,2,sum)>=3)] 

# load graph
network=HumanNet
Graph=graph.data.frame(network)

# find samples that have both mutation and expression data
com_samples=intersect(colnames(expr_data),mut_data_morethan_2)[1:10]
data_mut_com=mut_data[,com_samples]
data_expr_com=expr_data[,com_samples]

# calculate similarity matrix
correlation_samples=cor(data_expr_com)

# calculate co-mutated genes
co_mut=t(data_mut_com)%*%data_mut_com

# get importance score of the genes in each sample
Importance_Score=list()
for (i in 1:length(com_samples)){
  sample_i=com_samples[i]
  print(paste0(i," ",sample_i))
  importance_score=hygraph_randomwalk(sample_i,data_mut_com,co_mut,correlation_samples,Graph)
  Importance_Score[[i]]=importance_score

}
names(Importance_Score)=com_samples


}
