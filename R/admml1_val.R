#' admml1_val function
#'
#' This function allows to compute admm objective
#' @param theta vector
#' @param x data
#' @param alpha ADMM parameter
#' @param u vector
#' @param G matrix
#' @param Dk matrix
#' @param rho scalar
#' @export
#' admml1_val
#'
admml1_val = function(theta, x, alpha, u, G, Dk, rho) {
  arg = alpha - Dk %*% theta + u
  -logmltheta(theta, x, G) + (rho/2)*sum(arg^2)
}
