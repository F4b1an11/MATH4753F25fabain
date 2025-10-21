#' mychisim
#'
#' @param n1 sample size for each simulated sample
#' @param sigma1 true population standard deviation sigma used to form the chi square
#' @param mean1 true population mean used to generate data
#' @param iter number of monte carlo repitions
#' @param ymax upper limit of y axis for the histo
#' @param legend_pos takes in pos such as "topright", "bottomleft", "center",...
#' @param ...
#'
#' @returns A list with components:
#' \itemize{
#'   \item `w` – numeric vector of simulated chi-square statistics (length `iter`).
#'   \item `summary` – summary statistics for `w`.
#'   \item `sd` – standard deviation of `w`.
#'   \item `fun` – character label `"Chi-sq"`.
#' }
#' @export
#'
#' @examples
#' mychisim(iter=10000,ymax=0.15)
mychisim<-function(n1=10,sigma1=3,mean1=5,iter=1000,ymax=0.1,legend_pos = "topright", ...){    # adjust ymax to make graph fit
  y1=rnorm(n1*iter,mean=mean1,sd=sigma1)# generate iter samples of size n1

  data1.mat=matrix(y1,nrow=n1,ncol=iter,byrow=TRUE) # Each column is a sample size n1

  ssq1=apply(data1.mat,2,var) # ssq1 is s squared

  w=(n1-1)*ssq1/sigma1^2      #chi-sq stat

  hist(w,freq=FALSE, ylim=c(0,ymax), # Histogram with annotation
       main=substitute(paste("Sample size = ",n[1]," = ",n1," statistic = ",chi^2)),
       xlab=expression(paste(chi^2, "Statistic",sep=" ")), las=1)
  lines(density(w),col="Blue",lwd=3) # add a density plot
  curve(dchisq(x,n1-1),add=TRUE,col="Red",lty=2,lwd=3) # add a theoretical curve
  title=expression(chi^2==frac((n[1]-1)*s^2,sigma^2)) #mathematical annotation -see ?plotmath
  legend(legend_pos,c("Simulated","Theoretical"),col=c("Blue","Red"),lwd=4,lty=1:2,bty="n",title=title) # Legend #
  return(list(w=w,summary=summary(w),sd=sd(w),fun="Chi-sq")) # some output to use if needed
}

