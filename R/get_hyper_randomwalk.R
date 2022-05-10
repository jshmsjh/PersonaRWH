#' get_hyper_randomwalk
#' 
#' get_hyper_randomwalk uses the transition possibility matrix, the mutation data to perform random walk with restart on hypergraph
#' theta is the restart probability at every step of the random walk.
#' output of get_hyper_randomwalk is a list which combines the Importance of the genes and the distance of vectors between one iteration
#' 
#' @param P the transition possibility matrix
#' @param H adjacent matrix for i-th sample
#' @param theta the restart probability at every step of the random walk
#' @examples
#' @export get_hyper_randomwalk
get_hyper_randomwalk <-
  function(P,H,theta=0.8) {
    
    v0=rep(1/nrow(H),nrow(H))    #initial vector
    teleport=rep(1/nrow(H),nrow(H)) #
    
    Distance=c()
    vi=v0
    for(k in 1:10){
      if(k==1) print("1st iteration")
      else if(k==2) print("2nd iteration")
      else print(paste0(k,"th iteration"))
      vj=vi
      vi=theta*t(P)%*%vi+(1-theta)*teleport
      dis=sum(abs(vj-vi))
      Distance=append(Distance,dis)
    }
    Importance=data.frame(vi[order(vi,decreasing = T),])
    colnames(Importance)=c("Importance")
    return(Importance)
  }
