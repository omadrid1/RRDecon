#' admml1_grad  function
#'
#' This function allows to compute admm gradient
#' @param theta vector
#' @param x data
#' @param alpha ADMM parameter
#' @param u vector
#' @param G matrix
#' @param Dk matrix
#' @param rho scalar
#' @export
#' admml1_grad
#'
admml1_grad = function(theta, x, alpha, u, G, Dk, rho) {
  arg = Dk %*% theta - (alpha + u)
  -logmlthetagrad(theta, x, G) + rho*crossprod(Dk, arg)
}
