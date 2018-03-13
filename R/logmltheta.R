#' logmltheta function
#'
#' logarithm of log-likehood
#' @param theta vector
#' @param counts vector
#' @param G matrix
#' @export
#' logmltheta
#'
logmltheta = function(theta, counts, G)
{
  lambdamu = exp(theta);
  d = length(counts)

  lambday =  drop(G%*%lambdamu)
  loglike = sum(counts*log(lambday) - lambday)

  return(loglike)
}
