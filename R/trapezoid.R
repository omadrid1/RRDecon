#' trapezoid function
#'
#' trapezoind rule for integrals
#' @param grid  locations of points
#' @param f function values
#' @export
#' trapezoid
#'
trapezoid = function(grid, f)
{

  temp =  sum((grid[2:length(grid)] - grid[1:(length(grid)-1)])*(f[2:length(f)] + f[1:(length(f)-1)]   )/2)
  return(temp )
}
