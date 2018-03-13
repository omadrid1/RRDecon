#' dnormix function
#'
#' This function allows to evaluate the pdf of a mixture of normals
#' @param x locations where to evaluate the pdf
#' @param weights weights on each mixture component
#' @param means centers of mixture components
#' @param tau2 variances of mixture components
#' @export
#' dnormix
#'
dnormix  =  function(x,weights,mu,tau2)
{
  weights = weights/sum(weights)
  temp = rep(0,length(x))
  for(i in 1:length(mu))
  {
    temp = temp +  weights[i]*dnorm(x,mu[i],sqrt(tau2[i]))
  }
  return(temp)
}
