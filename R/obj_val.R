#' obj_val function
#'
#' objective function in L2 deconvolution
#' @param theta vector
#' @param x vector
#' @param G matrix
#' @param L matrix
#' @param lambda scalar
#' @export
#' obj_val
#'
obj_val = function(theta, x, G, L, lambda) {
  -logmltheta(theta, x, G) + (lambda/2)*{sum(theta * drop(L %*% theta))}
}
