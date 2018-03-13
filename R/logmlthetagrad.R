#' logmlthetagrad function
#'
#' logarithm of log-likehood
#' @param theta vector
#' @param counts vector
#' @param G matrix
#' @export
#' logmlthetagrad
#'
logmlthetagrad = function(theta,  counts,  G){
  d = length(counts)
  lmabday = rep(0,d)
  lambdamu = exp(theta);

  lambday = drop(G %*%  lambdamu)

  grad =  rep(0,d)
  grad =  drop(t(G) %*% (counts/lambday  - 1 ))
  grad = grad*  lambdamu


  return(grad)
}
