
#' @title Generate simulated data
#' @description The complete procedure for generating simulated data. For details, refer to simulation study (Section 3.1 in the main text).
#' @importFrom igraph as_adjacency_matrix sample_pa sample_gnp
#' @importFrom  MASS mvrnorm
#' @importFrom  stats rbeta rbinom rnorm runif var
#' @usage generate.data(p, n, K, network.type, umin, umax)
#' @param p The number of genes.
#' @param n The sample size. A positive intege or a vector with length equal to the number of cancer subtypes.
#' @param K The number of cancer subtypes.
#' @param network.type A character string indicating which network type is generated. "ER" (Erdos-Renyi) and "SF" (scale-free) can be used.
#' @param umin  The lower limits of the edge values.
#' @param umax  The upper limits of the edge values.
#' @details The function is used to generate the gene expression datasets.
#' @return
#' \item{\code{X}}{A list (length = \eqn{K}) of data matrices (\eqn{n_k \times p}), where K is the number of cancer subtypes and \eqn{n_k} is the sample size of k-th cancer subtype.}
#' \item{\code{theta.y}}{ A list (length = \eqn{K}) of precision matrices (\eqn{p \times p}) from cancerous cells of different  cancer subtyes, where K is the number of cancer subtypes
#'                       and \eqn{n_k} is the sample size of of k-th cancer subtype.}
#' \item{\code{theta.z}}{ A matrix of precision matrice  (\eqn{p \times p}) from non-cancerous cells shared by all cancer subtyes.}
#' \item{\code{purity}}{A list (length = \eqn{K}) of the tumor purity information vectors (\eqn{n_k \times 1}), where K is the number of cancer subtypes
#' and \eqn{n_k} is the sample size of of k-th cancer subtype.}
#'
#' @references  Jia-Juan Tu, Le Ou-Yang, Hong Yan, Xiao-Fei Zhang and Hong Qin (2019), Joint reconstruction of multiple gene networks by simultaneously capturing inter-tumor
#'      and intra-tumor heterogeneity.
#' @author Jia-Juan Tu
#' @seealso { \code{\link{NETI2}}, \code{\link{TCGA.BRCA}}}
#' @export
#' @examples
#' # Simulation data
#' data.x = generate.data(p = 100, n = 100, K = 4, network.type = "ER", umin = 0.5, umax = 1)


generate.data <- function(p = 100, n = 100, K = 4, network.type = "ER", umin = 0.5, umax = 1 ){

  if(length(n)==1){
    nvec=matrix(1,K,1)*n
  } else {
    nvec=n
  }

  Uy = matrix(list(), K,1)
  theta.y = matrix(list(), K,1)
  realP = matrix(list(), K,1)
  Y = matrix(list(), K,1)
  Z = matrix(list(), K,1)
  X = matrix(list(), K,1)
  purity = matrix(list(), K,1)

  if (network.type == "ER"){
    net  = igraph::as_adjacency_matrix(igraph::sample_gnp(p, 0.05), type = "both", sparse=FALSE)
  }

  if (network.type == "SF"){
    net = igraph::as_adjacency_matrix(igraph::sample_pa(p,directed = FALSE), type = "both", sparse=FALSE)
  }


  for (k in 1:K){
    temp = matrix(runif(p*p, umin, umax)*(2*rbinom(p*p, 1, 0.5) - 1), p, p)
    N = net*temp
    N = N*upper.tri(N)
    N = N + t(N)
    eigval_min = min(eigen(N)$values)
    theta.y[[k,1]] = N+ diag(p)*(abs(eigval_min) + 0.1)
    Uy[[k,1]] = rnorm(p,0,1)

  }


  temp = matrix(runif(p*p, umin, umax)*(2*rbinom(p*p, 1, 0.5) - 1), p, p)
  N = net*temp
  N = N*upper.tri(N)
  N = N + t(N)
  eigval_min = min(eigen(N)$values)
  theta.z = N+ diag(p)*(abs(eigval_min) + 0.1)
  Uz = rnorm(p,0,1)


  for (k in 1:K){
    Y[[k,1]] =  MASS::mvrnorm(nvec[k],Uy[[k,1]],solve(theta.y[[k,1]]))
    Z[[k,1]] =  MASS::mvrnorm(nvec[k],Uz,solve(theta.z))
  }


  meanP = 0.6;
  varP = 0.04;
  sizeP = meanP*(1-meanP)/varP-1
  for (k in 1:K){
    realP[[k,1]] = as.matrix(rbeta(nvec[k],meanP*sizeP,(1-meanP)*sizeP),ncol=1)


  }

  for (k in 1:K){
    X[[k,1]] = sweep(Y[[k,1]],1,realP[[k,1]],"*") + sweep(Z[[k,1]],1,1-realP[[k,1]],"*")
  }


  return(list(X = X,theta.y = theta.y, theta.z = theta.z, purity =  realP))
}

