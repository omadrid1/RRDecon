#' admml2_grad  function
#'
#' This function allows to compute admm gradient for L2 deconvol
#' @param theta vector
#' @param x data
#' @param alpha ADMM parameter
#' @param u vector
#' @param G matrix
#' @param Dk matrix
#' @param rho scalar
#' @export
#' admml2_grad
#'
admml2_grad = function(theta, x, b, G, D, L, alpha) {
  -logmlthetagrad(theta, x, G) + alpha*(L %*% theta - crossprod(D, b))
}
