#' admml2_val  function
#'
#' This function allows to compute admm objective for L2 deconvol
#' @param theta vector
#' @param x data
#' @param alpha ADMM parameter
#' @param u vector
#' @param G matrix
#' @param Dk matrix
#' @param rho scalar
#' @export
#' admml2_val
#'
admml2_val = function(theta, x, b, G, D, L, alpha) {
  arg = D %*% theta - b
  -logmltheta(theta, x, G) + (alpha/2)*sum(arg^2)
}
