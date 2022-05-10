\name{HumanNet}
\alias{HumanNet}
\docType{data}
\title{
HumanNet
}
\description{
HumanNet is a probabilistic functional gene network of 18,714 validated protein-encoding genes of Homo sapiens (by NCBI March 2007), constructed by a modified Bayesian integration of 21 types of 'omics' data from multiple organisms, with each data type weighted according to how well it links genes that are known to function together in H. sapiens. Each interaction in HumanNet has an associated log-likelihood score (LLS) that measures the probability of an interaction representing a true functional linkage between two genes.
}
\usage{data(HumanNet)}
\details{
Each interaction in HumanNet has an associated log-likelihood score (LLS) that measures the probability of an interaction representing a true functional linkage between two genes.
}
\source{
http://www.functionalnet.org/humannet/about.html
}
\references{
Insuk Lee, U. Martin Blom, Peggy I. Wang, Jung Eun Shin, and Edward M. Marcotte
Genome Research 21(7):1109-21 (2011)
}
\examples{
data(HumanNet)
}
\keyword{datasets}
