#' L2  deconconvolution function
#'
#' This function allows to perform L2 deconvolution with model selection as in the paper "A deconvolution path for mixtures"
#' @param y : a vector containing the raw observations.
#' @param prop : a scalar between 0 and 1  indicating the proportion of samples to be used as held out set. The default choice is 0.25
#' @param d : number of bins,  the default choice is floor(  (length(y)^(1/(2.01)) ))
#' @param lambda_grid : list of reuglarization parameters, default choice is  sort(c(0.001,0.1,5,seq(10,400,length = 16),seq(500,2000,length = 10),20000))
#' @return loc : location of the centers of bins where the density is estimated
#' @return f_hat : estimated mixing density
#' @return y_hat : marginal density
#' @export
#' L2_deconvolution
#'
L2_deconvolution = function(y,prop  =  0.25,d = floor(  (length(y)^(1/(2.01)) )) , lambda_grid = sort(c(0.001,0.1,5,seq(10,400,length = 16),seq(500,2000,length = 10),20000)))
{

  N = 1;
  #  prop = .35
  test_length = floor(prop*length(y))

  ycounts = rep(0, d)
  test_counts = rep(0,d)

  for(i in 1:N)
  {
    indexes =  sample.int(n, size = test_length , replace = FALSE, prob = NULL)
    p1 = gridprep(c(y[ -indexes ],min(y),max(y)),d)
    ycounts = ycounts + p1$ycounts/N # training counts

    test = y[indexes]
    test_hist = hist(test, p1$breaks,plot=FALSE)
    test_counts = test_counts + test_hist$counts/N
  }

  delta = mean(diff(p1$breaks))
  ycounts = floor(ycounts)
  test_counts = floor(test_counts)

  # Blur and penalty matrices
  G = as.matrix(p1$G)
  D = getDtf(d, 1)  # penalty matrix
  L = crossprod(D)

  # Initialization
  theta_hat = rep(log(mean(ycounts)), length(ycounts))

  # Fit via BFGS
  #sort(c(0.1,5,seq(10,400,length = 16),seq(500,2000,length = 10)))


  minus_test_likeilhood = 10^5

  error = rep(0,length(lambda_grid))
  cv = rep(0,length(error))

  for(k in 1:length(lambda_grid))
  {
    lambda = lambda_grid[k]

    system.time(out <- optim(theta_hat, obj_val, obj_grad, x = ycounts, G = G, L = L, lambda=lambda,
                             method='BFGS', control=list(maxit=10000, reltol=1e-12)))
    theta_hat = out$par
    f_hat = exp(theta_hat - log(n*delta))

    y_hat = G %*% f_hat

    ##########


    #     dof =  lambda*sum(abs(D%*% theta_hat) > 1e-2)
    theta_hat_test =   theta_hat - log(n*delta) + log(test_length*delta)

    zeta = D %*% theta_hat_test
    cv[k] =   -logmltheta(theta_hat_test, test_counts, G) + sum(abs(zeta)) #sum(abs(zeta) > 5e-4)


  }
  ind = which(cv == min(cv))




  best_lambda =  lambda_grid[ind]
  lambda =  best_lambda
  ########################
  p1 = gridprep(y,d)
  G = as.matrix(p1$G)
  ycounts = p1$ycounts
  lambda = best_lambda

  #     theta_hat = rep(log(mean(ycounts)), length(ycounts))
  #     f_true = dnormix(p1$mids,true_weights,true_mu,true_tau2)
  system.time(out <- optim(theta_hat, obj_val, obj_grad, x = ycounts, G = G, L = L, lambda=lambda,
                           method='BFGS', control=list(maxit=10000, reltol=1e-12)))
  theta_hat = out$par
  f_estimated = exp(theta_hat - log(n*delta))

  #   y_hat = G %*% f_estimated
  y_hat = f_estimated

  for(j in 1:d) {
    joint_dens = dnorm(p1$mids , p1$mids[j], 1) * f_estimated;
    y_hat[j] = trapezoid(p1$mids, joint_dens);
  }

  return(list(loc =  p1$mids, f_hat =f_estimated ,y_hat = y_hat ))
}

