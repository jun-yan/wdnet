#' Degree centrality
#'
#' Compute the degree centrality measure for a given weighted
#' directed network represented by its adjacency matrix. This is an
#' alternative to "degree_w" funciton in "tnet" package. Function
#' "degree_c" is matrix friendly, but "degree_w" function in "tnet"
#' package requires weight lists.
#'
#' @usage
#' degree_c(adj, alpha = 1, type = "out")
#'
#' @param adj Adjacency matrix of an weighted directed network
#' @param alpha tuning parameter; possible values
#' @param type which type to compute: "out" or "in"
#'
#' @return a scalar measure of centrality
#'
#' @references
#' Zhang, P. and Yan, J. (2020+) A review of centrality measure for
#'     weighted directed networks
#'
#' @examples
#' ## need to define adj_test first for the check to work
#' system.time(mydegree <- degree_c(adj_test, alpha = 0.8, type = "in"))
#' 
#' @export

degree_c <- function(adj, alpha = 1, type = "out") {
  if (alpha < 0){
    stop("The tuning parameter alpha must be nonnegative!")
  }
  if (dim(adj)[1]!=dim(adj)[2]) {
    stop("The adjacency matrix must be a square matrix!")
  }
  else{
    deg_c_output <- matrix(NA, nrow = dim(adj)[1], ncol = 2)
    deg_c_output[,1] <- c(1:dim(adj)[1])
    colnames(deg_c_output) <- c("vertex","degree_centrality")
    adj_deg <- adj
    adj_deg[which(adj_deg > 0)] <- 1
    if (type == "in"){
      deg_c_output[,2] <- rowSums(adj)^alpha + rowSums(adj_deg)^(1 - alpha)
    }
    if (type == "out"){
      deg_c_output[,2] <- colSums(adj)^alpha + colSums(adj_deg)^(1 - alpha)
    }
    return(deg_c_output)
  }
}

#' Closeness Centrality
#'
#' Compute closeness centrality values of vertices in a network by
#' introducing a tuning parameter. It is a matrix-friendly verstion of
#' "closeness_w" function in "tnet" package; the corresponding method
#' is "standard". The function allows for standard and harmonic
#' computations. It is equivalent to "closeness" function in "igraph"
#' package when method = "standard" and "alpha = 1"
#'
#' @param adj Adjacency matrix of a network
#' @param alpha numeric tuning parameter
#' @param type something
#' @param method something
#'
#' @return
#'
#' @examples
#' system.time(mydegree <- degree_c(adj_test, alpha = 0.8, type = "in"))
#' 
#' @export

closeness_c <- function(adj, alpha = 1, type = "out",
                        method = "harmonic") {
  ## require(igraph) # which function is from igraph? Use
  ## igraph::function
    
  if (alpha < 0) {
    stop("The tuning parameter alpha must be nonnegative!")
  }
  if (dim(adj)[1] != dim(adj)[2]) {
    stop("The adjacency matrix must be a square matrix!")
  }
  else {
    temp_g <- graph_from_adjacency_matrix(adj)
    closeness_c_output <- matrix(NA, nrow = dim(adj)[1], ncol = 2)
    closeness_c_output[,1] <- c(1:dim(adj)[1])
    colnames(closeness_c_output) <- c("vertex","closeness")
    if (method == "harmonic"){
      temp_d <- 1/distances(temp_g, mode = type, algorithm = "dijkstra")
      temp_d <- temp_d^alpha
      temp_d[temp_d == Inf] <- 0
      if (type == "in"){
        closeness_c_output[,2] <- colSums(temp_d)
      }
      if (type == "out"){
        closeness_c_output[,2] <- rowSums(temp_d)
      }
    }
    if (method == "standard"){
      temp_d <- distances(temp_g, mode = type, algorithm = "dijkstra")
      temp_d <- temp_d^alpha
      temp_d[temp_d == Inf] <- 0
      if (type == "in"){
        closeness_c_output[,2] <- 1/colSums(temp_d)
      }
      if (type == "out"){
        closeness_c_output[,2] <- 1/rowSums(temp_d)
      }
    }
    return(closeness_c_output)
    ## options(warn) = -1 # wrong place
  }
}


### Function "betweenness_c" computes betweenness centrality values of vertices in a network by introducing a tuning parameter
### Again, "betweenness_w" is a matrix-friendly verstion of "closeness_w" function in "tnet" package; the corresponding method is "standard"
### By default, network is set to be "directed," meaning that i -> j and j -> i are treated separately as two edges even in undirected networks. 
### The function is equivalent to "betweenness" function in "igraph" package when directed = "TRUE" and "alpha = 1"


betweenness_c <- function(adj, alpha = 1, directed = TRUE){
  require(tnet)
  if (alpha < 0){
    stop("The tuning parameter alpha must be nonnegative!")
  }
  if (dim(adj)[1]!=dim(adj)[2]){
    stop("The adjacency matrix must be a square matrix!")
  }
  if ((isSymmetric(adj) == FALSE) & (directed == FALSE)){
    stop("The adjacency matrix is not symmetric!")
  }
  if ((isSymmetric(adj) == TRUE) & (directed == TRUE)){
    warning("The adjacency matrix is symmetric!")
  }
  diag(adj) <- 0
  d <- dim(adj)[1]
  v.from <- rep(c(1:d),d)
  v.to <- sort(rep(c(1:d),d))
  v.weight <- as.vector(adj)
  temp_edgelist <- cbind(v.from, v.to, v.weight)
  if (directed == TRUE){
    suppressWarnings(tnet_edgelist <- as.tnet(temp_edgelist, type = "weighted one-mode tnet")) 
    myres <- betweenness_w(tnet_edgelist, directed = TRUE, alpha)
    colnames(myres) <- c("vertex","betweenness")
    return(myres)
  }
  if (directed == FALSE){
    suppressWarnings(tnet_edgelist <- as.tnet(temp_edgelist, type = "weighted one-mode tnet")) 
    myres <- betweenness_w(tnet_edgelist, directed = FALSE, alpha)
    colnames(myres) <- c("vertex","betweenness")
    return(myres)
  }
}

start_time_betweenness_c <- Sys.time()

mybetweenness <- betweenness_c(adj_test, alpha = 0.8, directed = TRUE)

end_time_betweenness_c <- Sys.time()

end_time_betweenness_c - start_time_betweenness_c


###--------------------------------------------------------------------------

###################
##### eigen_c #####
###################

### Function "eigen_c" computes eigenvector centrality values of vertices in a network
### It is equivalent to "eigen_centrality" function in "igraph" package, but it is much faster
### Besides, our function is matrix friendly, no need to convert to a graphic object
### The "eigen_c" function is based on "rARPACK" package maintained by Yixuan Qiu

eigen_c <- function(adj){
  require(rARPACK)
  if (dim(adj)[1]!=dim(adj)[2]){
    stop("The adjacency matrix must be a square matrix!")
  }
  if (dim(adj)[1] == 2){
    temp <- eigen(t(adj))
    eigen_v <- temp$vectors[,1]
    eigen_vstd <- abs(eigen_v)/max(abs(eigen_v))
    name_v <- c(1:dim(adj)[1])
    myres <- cbind(name_v, eigen_vstd)
    colnames(myres) <- c("vertex","eigenvec_centrality")
    return(myres)
  }
  if (dim(adj)[1] > 2){
    temp <-eigs(t(adj), k = 1, which = "LM", mattype = "matrix")
    eigen_v <- Re(temp$vectors)
    eigen_vstd <- abs(eigen_v) / max(abs(eigen_v))
    name_v <- c(1:dim(adj)[1])
    myres <- cbind(name_v, eigen_vstd)
    colnames(myres) <- c("vertex","eigenvec_centrality")
    return(myres)
  }
}

start_time_eigen_c <- Sys.time()

myeigen <- eigen_c(adj_test)

end_time_eigen_c <- Sys.time()

end_time_eigen_c - start_time_eigen_c

###--------------------------------------------------------------------------

##################
##### katz_c #####
##################

### Function "katz_c" computes the Katz centrality values of vertices in a network
### It is equivalent to "alpha_centrality" function in "igraph" package, but it is more general
### The function "alpha_centrality" only allows constant exogenous values for all the vertices
### Our algorithm allows for a more general class

katz_c <- function(adj, alpha, beta = rep(1, dim(adj)[1])){
  require(rARPACK)
  if (dim(adj)[1]!=dim(adj)[2]){
    stop("The adjacency matrix must be a square matrix!")
  }
  if (length(beta)!=dim(adj)[1]){
    stop("The dimensions of beta and the adjacency matrix are not equal!")
  }
  n <- dim(adj)[1]
  if (n == 2){
    spec <- 1/Re(eigen(adj)$value)  
  } else {
    spec <- 1/Re(eigs(adj, k = 1, which = "LM", mattype = "matrix")$values)
  }
  if (alpha > spec){
    warning("Alpha is greater than the spectral radius of the adjacency matrix!")
  }
  temp <- diag(n) - alpha*t(adj)
  temp_inv <- solve(temp)
  katz_v <- temp_inv %*% beta
  name_v <- c(1:dim(adj)[1])
  myres <- cbind(name_v, katz_v)
  colnames(myres) <- c("vertex","Katz_centrality")
  return(myres)
}

start_time_katz_c <- Sys.time()

mykatz <- katz_c(adj_test, alpha = 4.5*10^(-5))

end_time_katz_c <- Sys.time()

end_time_katz_c - start_time_katz_c

###--------------------------------------------------------------------------

####################
##### PageRank #####
####################

### Function "pagerank_c" computes PageRank centrality values of vertices in a network
### It is an extension of "page_rank" function in "igraph" package, but it is much faster
### Our function is matrix friendly, no need to convert to a graphic object
### Another extension is that our function allows for prior information if it is avaialble
### The "pagerank_c" function is based on "rARPACK" package maintained by Yixuan Qiu

pagerank_c <- function(adj, alpha = 0.85, prior.info = rep(1/dim(adj)[1],dim(adj)[1])){
  require(rARPACK)
  if ((alpha < 0) | (alpha > 1)){
    stop("The damping factor must be between 0 and 1!")
  }
  if (dim(adj)[1]!=dim(adj)[2]){
    stop("The adjacency matrix must be a square matrix!")
  }
  if (length(prior.info)!=dim(adj)[1]){
    stop("The dimensions of the prior information and the adjacency matrix are not equal!")
  }
  if (dim(adj)[1] == 2){
    sinks <- which(rowSums(adj) == 0)
    M <- matrix(NA, 2, 2)
    M <- adj / rowSums(adj)
    M[sinks,] <- prior.info
    E <- matrix(1, 2, 2)
    Mstar <- alpha*t(M) + (1 - alpha)/2*E 
    temp <- eigen(Mstar)
    eigen_v <- temp$vectors[,1]
    eigen_vstd <- abs(eigen_v)/sum(abs(eigen_v))
    name_v <- c(1:dim(adj)[1])
    myres <- cbind(name_v, eigen_vstd)
    colnames(myres) <- c("vertex","PageRank")
    return(myres)
  }
  if (dim(adj)[1] > 2){
    n <- dim(adj)[1]
    sinks <- which(rowSums(adj) == 0)
    M <- matrix(NA, n, n)
    M <- adj / rowSums(adj)
    M[sinks,] <- prior.info
    E <- matrix(1, n, n)
    Mstar <- alpha*t(M) + (1 - alpha)/n*E 
    temp <-eigs(Mstar, k = 1, which = "LM", mattype = "matrix")
    eigen_v <- Re(temp$vectors)
    eigen_vstd <- abs(eigen_v) / sum(abs(eigen_v))
    name_v <- c(1:dim(adj)[1])
    myres <- cbind(name_v, eigen_vstd)
    colnames(myres) <- c("vertex","PageRank")
    return(myres)
  }
}

start_time_pagerank_c <- Sys.time()

mypagerank <- pagerank_c(adj_test)

end_time_pagerank_c <- Sys.time()

end_time_pagerank_c - start_time_pagerank_c

###--------------------------------------------------------------------------

################
##### HITS #####
################

### Function "HITS_c" computes hub scores and authorities scores of vertices in a network
### This is a matrix friendly version of "hub_score" and "authority_score" functions in igraph package

HITS_c <- function(adj){
  require(rARPACK)
  if (dim(adj)[1]!=dim(adj)[2]){
    stop("The adjacency matrix must be a square matrix!")
  }
  
  A <- adj%*%t(adj)
  
  if (dim(A)[1] == 2){
    temp <- eigen(A)
    eigen_v <- temp$vectors[,1]
    eigen_vstd <- abs(eigen_v)/max(abs(eigen_v))
    name_v <- c(1:dim(A)[1])
    eigen_t <- t(adj)%*%eigen_vstd
    eigen_tstd <- abs(eigen_t) / max(abs(eigen_t))
    myres <- cbind(name_v, eigen_vstd,eigen_tstd)
    colnames(myres) <- c("vertex","hubs","authorities")
    return(myres)
  }
  if (dim(A)[1] > 2){
    temp <-eigs(t(A), k = 1, which = "LM", mattype = "matrix")
    eigen_v <- Re(temp$vectors)
    eigen_vstd <- abs(eigen_v) / max(abs(eigen_v))
    name_v <- c(1:dim(A)[1])
    eigen_t <- t(adj)%*%eigen_vstd
    eigen_tstd <- abs(eigen_t) / max(abs(eigen_t))
    myres <- cbind(name_v, eigen_vstd,eigen_tstd)
    colnames(myres) <- c("vertex","hubs","authorities")
    return(myres)
  }
}

start_time_eigen_c <- Sys.time()

myeigen <- eigen_c(adj_test)

end_time_eigen_c <- Sys.time()

end_time_eigen_c - start_time_eigen_c

