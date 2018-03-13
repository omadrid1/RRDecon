#' soft function
#'
#' soft thresholding
#' @param x vector
#' @param lambda penalty
#' @export
#' soft
#'
soft = function(x, lambda) {
  sign(x) * pmax(0, abs(x) - lambda)
}
