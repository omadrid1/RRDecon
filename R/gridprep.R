#' gridprep function
#'
#' Preparing data for calculations on grid
#' @param y raw observations for deconvolution
#' @param gridsize  this is the number of bins
#' @export
#' gridprep
#'
gridprep = function(y, gridsize=250) {
  n = length(y)

  # Set up the discrete grid
  ymin = min(y) - 1
  ymax = max(y) + 1
  mybreaks = seq(ymin, ymax, length=gridsize+1)
  k = length(mybreaks)
  mids = (mybreaks[2:k] + mybreaks[1:(k-1)])/2

  yhist = hist(y, mybreaks,plot=FALSE)
  ycounts = yhist$counts
  G = outer(mids,mids, FUN=function(x,y){dnorm(x-y)*diff(mybreaks)})

  list(mids=mids, ycounts=ycounts, breaks=mybreaks, G=G)
}
