#' ntickets
#'
#' @param N # of seats
#' @param gamma probability of overbooking
#' @param p probability of passenger showing up
#'
#' @returns list(nd,nc, N,p,gamma) nd is the optimal number of seats sold (discrete).
#' nc is the optimal number of seats sold (continuous)
#' @export
#'
#' @examples
#' ntickets(200,0.02,0.95)
ntickets<-function(N, gamma,p){
  # N = # of Seats
  # n = # of tickets sold
  # y = probability of overbooking
  # p = probability of passenger showing up

  # --- Discrete (exact) calculation ---
  # Find smallest n such that P(X > N) <= gamma, X ~ Binomial(n, p)
  xx<- seq(N,N*1.1, by = 1)
  obj <- pbinom(N,xx,p)-1 + gamma
  index <- which.min(abs(obj))
  nd <- xx[index]

  # --- Normal approximation calculation ---
  # N + 0.5 = qnorm(1 - gamma, np, sqrt(np(1-p)))
  # N + 0.5 = np + z*sqrt(np(1-p)), z = omega^-1(1-gamma)
  n_from_qnorm <- function(N,gamma,p){
  z <- qnorm(1-gamma)
  a <- p
  b <- sqrt(p*(1-p))
  t <- (-z*b+ sqrt((z*b)^2 + 4*a*(N+0.5)))/(2*a)
  nc <- t^2
  return(nc)
  }
  nc <- n_from_qnorm(N,gamma,p)

  plot(xx,abs(obj),
       pch = 21,
       xlab = "Number of Tickets sold, n",
       ylab = "Size of the objective function",
       bg = ifelse(xx != nd, "blue", "red"),
       type = "b",
       main = paste("Objective Vs n to find optimal tickets sold (", nd,
         ")\n gamma = ", gamma, " N = ", N, " discrete")
       )
  abline(h = abs(pbinom(N,nd,p)-1 + gamma),v = nd, col = "red")

  curve(xlim = c(N,N*1.1),
        abs(pnorm((N+0.5 - x*p)/sqrt(x*p*(1-p)))-1 + gamma),
        xlab = "Number of Tickets sold, n",
        ylab = "Size of the objective function",
        main = paste("Objective Vs n to find optimal tickets sold (", nc,
                     ")\n gamma = ", gamma, " N = ", N, " continous")
        )
  abline(h = pnorm((N+0.5 - nc*p)/sqrt(nc*p*(1-p)))-1 + gamma, v = nc, col = "red")

  list(nd = nd ,nc = nc ,N = N,p = p,
       gamma = gamma
       )
}













