#' L1  deconconvolution function
#'
#' This function allows to perform L1 deconvolution with model selection as in the paper "A deconvolution path for mixtures"
#' @param y : a vector containing the raw observations.
#' @param d : number of bins
#' @return loc : location of the centers of bins where the density is estimated
#' @return f_hat : estimated mixing density
#' @return y_hat : marginal density
#' @export
#' L1_deconvolution
#'

L1_deconvolution =  function(y,d)
{
# source("logmlthetagrad.R")
 # source("logmltheta.R")
#  source("Dtf2_prod_x.R")

  # Set up the problem
  p1 = gridprep(y,d)
  ycounts = p1$ycounts
  delta = mean(diff(p1$breaks))

  # Blur and penalty matrices
  G = as.matrix(p1$G)
  K = 2
  D = getDtf(d, K) # penalty matrix
  L = crossprod(D)

  # Initialization
  theta_hat = rep(log(mean(ycounts)), length(ycounts))

  # Input pars
  #  d = length(ycounts)
  K = 1
  abs_tol = 5e-4
  inflate = 1.5
  n_lambda = 15

  # Matrices
  D1 = getDtf(d-K, 0) # penalty matrix
  Dk = getDtf(d, K-1)
  Dtf = D1 %*% Dk
  DtD = crossprod(Dk)

  # A grid of lambda values
  lambda_grid = 10^seq(3,-1,length=n_lambda)

  dev_trace = rep(0, n_lambda)
  dof_trace = rep(0, n_lambda)
  aic_trace = rep(0, n_lambda)
  theta_trace = matrix(0, nrow=d, ncol=n_lambda)

  # Initialize
  theta = rep(log(mean(ycounts)), length(ycounts))
  u = rep(0, nrow(Dk))
  lambda = lambda_grid[1]
  rho = lambda
  alpha = glmgen::trendfilter(drop(Dk %*% theta - u), k = 0L, lambda=lambda/rho)$beta[,1]

  # Solution path
  for(i in seq_along(lambda_grid)) {
    lambda = lambda_grid[i]
    #  cat('lambda =', lambda, '\n')
    rho = lambda
    step_counter = 0
    converged=FALSE
    while(!converged && step_counter < 1000) {
      step_counter = step_counter + 1

      ## Two-stage quasi Newton step
      # First, warm-start theta using BFGS
      out <- optim(theta, admml1_val, admml1_grad,
                   x = ycounts, alpha=alpha, u = u, G = G, Dk=Dk, rho=rho,
                   method='BFGS', control=list(maxit=10000, reltol = 1e-21))
      theta = out$par


      # Update alpha
      diff_theta = Dk %*% theta
      alpha_arg = drop(diff_theta - u)
      alpha_new = glmgen::trendfilter(alpha_arg, k = 0L, lambda=lambda/rho)$beta[,1]
      dual_residual = rho*(alpha_new - alpha)
      alpha = alpha_new

      # Update u
      primal_residual = drop(alpha - diff_theta)
      u = u + primal_residual

      # Check convergence
      primal_resnorm = sqrt(mean(primal_residual^2))
      dual_resnorm = sqrt(mean(dual_residual^2))
      if(dual_resnorm < abs_tol && primal_resnorm < abs_tol) {
        converged=TRUE
      }

    }

    zeta = Dtf %*% theta
    dev = -logmltheta(theta, ycounts, G)
    dof = sum(abs(zeta) > 5e-4)
    aic = 2*dev + 2*dof
    dev_trace[i] = dev
    dof_trace[i] = dof
    aic_trace[i] = aic
    theta_trace[,i] = theta
  }###  clsoing loop for solution path

  jbest = which.min(aic_trace)
  theta = theta_trace[,jbest]
  zeta = Dtf %*% theta

  theta_hat = theta
  # f_hat = exp(theta_hat - log(n*delta)

  #f_true = dnormix(p1$mids,true_weights,true_mu,true_tau2)

  f_estimated = exp(theta_hat - log(n*delta))
  y_hat = G %*%  f_estimated

  y_hat = f_estimated

  for(j in 1:d) {
    joint_dens = dnorm(p1$mids , p1$mids[j], 1) * f_estimated;
    y_hat[j] = trapezoid(p1$mids, joint_dens);
  }

  return(list(loc  =  p1$mids,f_hat =f_estimated ,y_hat = y_hat))
}###  clsoing function


