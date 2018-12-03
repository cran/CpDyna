
#' Function to simulate a dataset.
#' @param Lbd A (K x K x D) array with Poisson stepwise intensities, where K is the number of clusters and D is the number of time segments. 
#' @param Pi A vector of length K. It contains the proportions of nodes in each cluster.
#' @param eta  A vector of length D containing the positive ordered change points. The last entry is equal to T.
#' @param N  An integer. The number of nodes in the graph.
#' @return \item{events}{ A (3 x Nb.events) matrix. The first row contains the ordered interaction times, the second row the source nodes, the third row the destination nodes. Currently, only undirected interactions are taken into account. Thus, for each interaction time, the source node is stricly smaller than the destination node. No self loops allowed.}
#' @return \item{z}{An integer vector of length N. It contains the cluster labels of nodes.}
ModSim <- function(Lbd,  # A (K x D) array with Poisson stepwise intensities
                   Pi,   # proportion of nodes to each class
                   eta,  # D positive change points, the last one equal to T
                   N     # number of nodes
){
  D <- length(eta)
  K <- nrow(Lbd)
  z <- sample(K, N, replace = TRUE, prob = Pi)
  events <- c()
  for ( i in 1:(N-1) ){
    for (j in (i+1):N){
      tmp <- c()
      for(d in 1:D){
        if (d == 1){
          lbd <- Lbd[z[i],z[j],d]*eta[d]
          Mij <- rpois(1, lbd)
          tmp <- cbind(tmp, matrix(c(sort(runif(Mij, min = 0, max = eta[d])), rep(i, times = Mij), rep(j, times = Mij)), nrow = 3, byrow = TRUE))    
        }
        else{ 
          lbd <- Lbd[z[i], z[j], d]*(eta[d] - eta[d-1])
          Mij <- rpois(1, lbd)
          tmp <- cbind(tmp, matrix(c(sort(runif(Mij, min = eta[d-1], max = eta[d])), rep(i, times = Mij), rep(j, times = Mij)), nrow = 3, byrow = TRUE))    
        }
      }
      events <- cbind(events, tmp) 
    } 
  }
  
  M <- ncol(events)
  pos <- events[1,]
  names(pos) <- c(1:M)
  pos <- sort(pos)
  pos <- as.numeric(names(pos))
  events <- events[,pos]
  
  return(list(events = events, z = z))
}
