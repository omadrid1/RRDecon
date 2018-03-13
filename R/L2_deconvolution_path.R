#' L2  deconconvolution path function
#'
#' This function allows to compute the deconvolution path using L2  regularization
#' @param y : a vector containing the raw observations.
#' @param lambda_grid : a vector containing a list of regularization parameters to compute the path. Deafault choice is sort(c(0.001,0.1,5,seq(10,400,length = 16),seq(500,2000,length = 10),20000)).
#' @param d : number of bins, the default  choice is floor(  (length(y)^(1/(2.01)) ))
#' @return loc : location of the centers of bins where the density is estimated
#' @return f_hat : a matrix containing the solution path of mixing densities, each row corresponds to regularization parameter.
#' @return y_hat : a matrix containing the solution path of marginal densities, each row corresponds to regularization parameter.
#' @export
#' L2_deconvolution_path

L2_deconvolution_path = function(y,lambda_grid = sort(c(0.001,0.1,5,seq(10,400,length = 16),seq(500,2000,length = 10),20000)),d = floor(  (length(y)^(1/(2.01)) )) )
{

  N = 1;

  ycounts = rep(0, d)
  p1 = gridprep(y,d)
  ycounts =  p1$ycounts

  delta = mean(diff(p1$breaks))
  ycounts = floor(ycounts)

  # Blur and penalty matrices
  G = as.matrix(p1$G)
  D = getDtf(d, 1)  # penalty matrix
  L = crossprod(D)

  # Initialization
  theta_hat = rep(log(mean(ycounts)), length(ycounts))

  # Fit via BFGS

  minus_test_likeilhood = 10^5

  f_hat   =  matrix(0, length(lambda_grid),d)
  y_hat   =  matrix(0, length(lambda_grid),d)

  for(k in 1:length(lambda_grid))
  {
    lambda = lambda_grid[k]

    system.time(out <- optim(theta_hat, obj_val, obj_grad, x = ycounts, G = G, L = L, lambda=lambda,
                             method='BFGS', control=list(maxit=10000, reltol=1e-12)))
    theta_hat = out$par
    f_hat[k,] = exp(theta_hat - log(n*delta))

    y_hat[k,] = G %*% f_hat[k,]

    f_estimated = f_hat[k,]
    #y_hat = f_estimated

    for(j in 1:d) {
      joint_dens = dnorm(p1$mids , p1$mids[j], 1) * f_estimated;
      y_hat[k,j] = trapezoid(p1$mids, joint_dens);
    }

    ##########

  }
  ########################
#
#   y_hat = f_estimated
#
#   for(j in 1:d) {
#     joint_dens = dnorm(p1$mids , p1$mids[j], 1) * f_estimated;
#     y_hat[j] = trapezoid(p1$mids, joint_dens);
#   }

  return(list(loc =  p1$mids, f_hat =f_hat ,y_hat = y_hat ))
}

