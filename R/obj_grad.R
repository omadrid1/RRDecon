#' obj_grad  function
#'
#' compute gradient of objective in L2 deconvolution
#' @param theta vector
#' @param x vector
#' @param G matrix
#' @param L matrix
#' @param lambda scalar
#' @export
#' obj_grad
#'
obj_grad = function(theta, x, G, L, lambda) {
  -logmlthetagrad(theta, x, G) + lambda*(L %*% theta)
}
