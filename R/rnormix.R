#' rnormix function
#'
#' This function allows to generate sample from a mixture of normals
#' @param n number of samples
#' @param weights weights on each mixture component
#' @param means centers of mixture components
#' @param tau2 variances of mixture components
#' @export
#' rnormix
#'

rnormix =  function(n, weights, means, tau2)
{
  weights = weights/sum(weights)
  cc =  sample.int(length(means), size = n, replace = TRUE, prob = weights)

  return( rnorm(n,means[cc], sqrt(tau2[cc])) )
}
