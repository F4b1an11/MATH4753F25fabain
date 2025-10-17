#' Hypergeometric Simulation
#'
#' Simulates repeated draws from a hypergeometric distribution by
#' sampling without replacement from a finite population.
#'
#' @param iter Number of simulation iterations to run. Default is 100.
#' @param N Total population size (number of marbles in the bag).
#' @param r Number of "success" items in the population
#'   (e.g., number of white marbles).
#' @param n Sample size drawn without replacement from the population.
#'
#' @return A table of relative frequencies (proportions) for the number
#'   of successes observed across the \code{iter} simulations.
#'   A barplot of the distribution is also displayed.
#'
#' @export
#'
#' @examples
#' # Run 100 simulations of drawing 5 marbles from a bag with
#' # 12 white marbles (successes) and 8 black marbles (failures)
#' myhyper(iter = 100, N = 20, r = 12, n = 5)
#'
#' # Larger simulation with 1000 iterations and sample size 19
#' myhyper(iter = 1000, N = 20, r = 12, n = 19)
#'
#' # Compare to theoretical hypergeometric probabilities
#' dhyper(x = 0:19, m = 12, n = 8, k = 19)
myhyper=function(iter=100,N=20,r=12,n=5){
  # make a matrix to hold the samples
  #initially filled with NA's
  sam.mat=matrix(NA,nr=n,nc=iter, byrow=TRUE)
  #Make a vector to hold the number of successes over the trials
  succ=c()
  for( i in 1:iter){
    #Fill each column with a new sample
    sam.mat[,i]=sample(rep(c(1,0),c(r,N-r)),n,replace=FALSE)
    #Calculate a statistic from the sample (this case it is the sum)
    succ[i]=sum(sam.mat[,i])
  }
  #Make a table of successes
  succ.tab=table(factor(succ,levels=0:n))
  #Make a barplot of the proportions
  barplot(succ.tab/(iter), col=rainbow(n+1), main="HYPERGEOMETRIC simulation", xlab="Number of successes")
  succ.tab/iter
}
