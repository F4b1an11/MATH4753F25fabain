#' Maximum Likelihood Estimation via Grid Search
#'
#' \code{mymaxlik()} evaluates a log-likelihood function over a grid of parameter
#' values and returns the parameter that maximizes the likelihood. It produces a
#' plot of the log-likelihood curve and marks the maximum likelihood estimate
#' (MLE) visually. This function is useful when analytical maximization is not
#' possible and numerical methods such as grid approximation are required.
#'
#' @param lfun A function that accepts two arguments \code{x} and \code{param}
#'   and returns the log-likelihood contribution for a single observation. This
#'   function must be vectorized via \code{outer()} to handle all combinations of
#'   data and parameter values.
#'
#' @param x A numeric vector containing the observed data.
#'
#' @param param A numeric vector of parameter values over which the likelihood
#'   should be evaluated. This forms the grid for the maximization process.
#'
#' @param ... Additional graphical parameters passed to \code{plot()},
#'   such as \code{main}, \code{xlab}, or \code{ylab}.
#'
#' @returns A list with the following components:
#' \item{i}{The index in \code{param} corresponding to the maximum likelihood estimate.}
#' \item{parami}{The parameter value at which the likelihood is maximized (the MLE).}
#' \item{yi}{The maximum value of the log-likelihood.}
#' \item{slope}{A numeric vector describing approximate slopes around the MLE,
#'   useful for verifying sign changes that indicate a maximum.}
#'
#' @export
#'
#' @examples
#' # Example: MLE for p in a Binomial(10, p) model
#' logbin <- function(x, p) log(dbinom(x, size = 10, prob = p))
#' x <- c(9, 9, 1, 9, 9, 9)
#' param <- seq(0, 1, length = 1000)
#' mymaxlik(lfun = logbin, x = x, param = param,
#'          xlab = expression(p), main = "Binomial MLE")
#'
mymaxlik=function(lfun,x,param,...){
  # how many param values are there?
  np=length(param)
  # outer -- notice the order, x then param
  # this produces a matrix -- try outer(1:4,5:10,function(x,y) paste(x,y,sep=" "))   to understand
  z=outer(x,param,lfun)
  # z is a matrix where each x,param is replaced with the function evaluated at those values
  y=apply(z,2,sum)

  # y is a vector made up of the column sums
  # Each y is the log lik for a new parameter value
  plot(param,y,col="Blue",type="l",lwd=2,...)
  # which gives the index for the value of y == max.
  # there could be a max between two values of the parameter, therefore 2 indices
  # the first max will take the larger indice
  i=max(which(y==max(y)))
  abline(v=param[i],lwd=2,col="Red")

  # plots a nice point where the max lik is
  points(param[i],y[i],pch=19,cex=1.5,col="Black")
  axis(3,param[i],round(param[i],2))
  #check slopes. If it is a max the slope shoud change sign from + to
  # We should get three + and two -vs
  ifelse(i-3>=1 & i+2<=np, slope<-(y[(i-2):(i+2)]-y[(i-3):(i+1)])/(param[(i-2):(i+2)]-param[(i-3):(i+1)]),slope<-"NA")
  return(list(i=i,parami=param[i],yi=y[i],slope=slope))
}
