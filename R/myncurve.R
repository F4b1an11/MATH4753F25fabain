#' myncurve creates a curve of an norm distribution
#'
#' @param mu is the mean, a real number
#' @param sigma is the standard deviation, a real number
#' @param a for P(Y < a), a real number
#'
#' @returns a list of mu and sigma
#' @export
#'
#' @examples
#' myncurve(1,1,1)
myncurve = function(mu, sigma, a){
  curve(dnorm(x,mean=mu,sd=sigma), xlim = c(mu-3*sigma, mu + 3*sigma))
  xcurve <- seq(-mu-3*sigma, a, length = 1000)
  ycurve <- dnorm(xcurve,mean = mu, sd = sigma)
  polygon(x = c((mu-3*sigma),xcurve,a) , y = c(0,ycurve,0), col = "orange")
  text(x = mu, y = 0.5*dnorm((a-mu-3*sigma)/2,mu,sigma),
       label = paste("Area = ",round(pnorm(a,mu,sigma),4) ))
  list(mu = mu, sigma = sigma)
}
myncurve(1,1,1)
