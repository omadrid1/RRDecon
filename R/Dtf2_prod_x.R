#' Dtf2_prod_x
#'
#' just auxiliar function
#' @param x vector
#' @export
#' Dtf2_prod_x
#'
Dtf2_prod_x = function(x)
{
  n = length(x)
  result = rep(0,n-3)
  result =  x[1:(n-3)] - 3 *x[ 2 :(n-2)] + 3*x[ 3:(n-1)] - x[4:n]
  #for(int i=0; i<(n-3); i++) {
   #  result[i] = x[i] - 3.0*x[i+1] + 3.0*x[i+2] - x[i+3];
  #}
  return(result)
}
