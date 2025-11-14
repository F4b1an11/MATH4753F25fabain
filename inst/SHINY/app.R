# app.R
# Shiny app: Demonstrating Maximum Likelihood Estimation for 5+ univariate distributions
# Distributions: Normal, Exponential, Poisson, Binomial, Gamma
# Author: You :)

library(shiny)

# ---------- Helper: numerical stability ----------
neg_inf <- -1e20

# ---------- Log-likelihoods ----------
ll_norm <- function(par, x) {
  # par = c(mu, sigma_log) to enforce sigma > 0
  mu <- par[1]
  sigma <- exp(par[2])
  if (sigma <= 0) return(neg_inf)
  sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
}

ll_exp <- function(par, x) {
  # par = c(log_lambda) to enforce lambda > 0
  lambda <- exp(par[1])
  if (lambda <= 0) return(neg_inf)
  sum(dexp(x, rate = lambda, log = TRUE))
}

ll_pois <- function(par, x) {
  # par = c(log_lambda) to enforce lambda > 0
  lambda <- exp(par[1])
  if (lambda <= 0) return(neg_inf)
  sum(dpois(x, lambda = lambda, log = TRUE))
}

ll_binom <- function(par, x, m) {
  # par = c(logit_p) to enforce 0 < p < 1
  p <- 1/(1 + exp(-par[1]))
  if (p <= 0 || p >= 1) return(neg_inf)
  sum(dbinom(x, size = m, prob = p, log = TRUE))
}

ll_gamma <- function(par, x) {
  # par = c(log_k, log_theta) shape k > 0, scale theta > 0
  k <- exp(par[1]); theta <- exp(par[2])
  if (k <= 0 || theta <= 0) return(neg_inf)
  sum(dgamma(x, shape = k, scale = theta, log = TRUE))
}

# ---------- Closed-form MLEs + SEs where available ----------
mle_closed <- function(dist, x, m_trials = NULL) {
  n <- length(x)
  out <- list(est = NULL, se = NULL, ll = NA, method = "closed-form")
  if (dist == "Normal") {
    mu_hat <- mean(x)
    # MLE sigma is sqrt( (1/n) * sum (x - mean)^2 ) (note: not sample sd with n-1)
    sigma_hat <- sqrt(mean((x - mu_hat)^2))
    out$est <- c(mu = mu_hat, sigma = sigma_hat)
    # Approx SEs: SE(mu_hat) = sigma / sqrt(n); SE(sigma_hat) ~ sigma/sqrt(2n)
    out$se  <- c(mu = sigma_hat/sqrt(n), sigma = sigma_hat/sqrt(2*n))
    out$ll  <- ll_norm(c(mu_hat, log(sigma_hat)), x)
  } else if (dist == "Exponential") {
    lambda_hat <- 1/mean(x)
    out$est <- c(lambda = lambda_hat)
    # Fisher info: I(λ) = n/λ^2  => Var(λ_hat) ~ λ^2 / n
    out$se  <- c(lambda = lambda_hat / sqrt(n))
    out$ll  <- ll_exp(c(log(lambda_hat)), x)
  } else if (dist == "Poisson") {
    lambda_hat <- mean(x)
    # SE(lambda_hat) = sqrt(lambda / n)
    out$est <- c(lambda = lambda_hat)
    out$se  <- c(lambda = sqrt(lambda_hat / n))
    out$ll  <- ll_pois(c(log(lambda_hat)), x)
  } else if (dist == "Binomial") {
    stopifnot(!is.null(m_trials))
    p_hat <- mean(x)/m_trials
    p_hat <- min(max(p_hat, 1e-8), 1 - 1e-8) # guard
    # SE(p_hat) = sqrt( p(1-p) / (n*m) ) for sample mean of Bin(m,p)
    out$est <- c(p = p_hat)
    out$se  <- c(p = sqrt(p_hat*(1 - p_hat)/(n*m_trials)))
    out$ll  <- ll_binom(c(qlogis(p_hat)), x, m_trials)
  } else {
    stop("No closed form for this distribution here.")
  }
  out
}

# ---------- Generic optimizer wrapper ----------
mle_optim <- function(llfun, par0, x, transform = identity, extra = NULL) {
  fn <- function(p) {
    if (is.null(extra)) return(-llfun(p, x))
    do.call(llfun, c(list(p, x), extra))
    return(-do.call(llfun, c(list(p, x), extra)))
  }
  opt <- optim(par = par0, fn = function(p) {
    if (is.null(extra)) return(-llfun(p, x))
    return(-do.call(llfun, c(list(p, x), extra)))
  }, control = list(reltol = 1e-10), method = "BFGS")
  list(par = opt$par, ll = -opt$value, ok = (opt$convergence == 0))
}

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("MLE Demo for Univariate Distributions"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dist", "Distribution:",
                  c("Normal", "Exponential", "Poisson", "Binomial", "Gamma")),
      sliderInput("n", "Sample size n:", min = 10, max = 5000, value = 200, step = 10),
      numericInput("seed", "Random seed:", value = 1, min = 1),
      conditionalPanel(
        condition = "input.dist == 'Normal'",
        numericInput("norm_mu", "True μ:", value = 0),
        numericInput("norm_sigma", "True σ (>0):", value = 1, min = 0.0001, step = 0.1)
      ),
      conditionalPanel(
        condition = "input.dist == 'Exponential'",
        numericInput("exp_lambda", "True λ (>0):", value = 1, min = 0.0001, step = 0.1)
      ),
      conditionalPanel(
        condition = "input.dist == 'Poisson'",
        numericInput("pois_lambda", "True λ (>0):", value = 5, min = 0.0001, step = 0.1)
      ),
      conditionalPanel(
        condition = "input.dist == 'Binomial'",
        numericInput("binom_m", "Trials m (integer ≥1):", value = 10, min = 1, step = 1),
        numericInput("binom_p", "True p in (0,1):", value = 0.3, min = 0.0001, max = 0.9999, step = 0.01)
      ),
      conditionalPanel(
        condition = "input.dist == 'Gamma'",
        numericInput("gamma_k", "True shape k (>0):", value = 2, min = 0.1, step = 0.1),
        numericInput("gamma_theta", "True scale θ (>0):", value = 2, min = 0.1, step = 0.1)
      ),
      actionButton("go", "Generate & Fit MLE", class = "btn-primary")
    ),
    mainPanel(
      h4("True parameters vs MLE"),
      tableOutput("est_table"),
      h4("Log-likelihood at MLE (and at True parameters)"),
      verbatimTextOutput("ll_text"),
      h4("Data & Fits"),
      plotOutput("fitPlot", height = "360px"),
      conditionalPanel(
        condition = "input.dist == 'Normal' || input.dist == 'Gamma'",
        h4("Log-likelihood contour (2-parameter family)"),
        plotOutput("llContour", height = "360px")
      ),
      hr(),
      p(em("Notes:")),
      tags$ul(
        tags$li("For Normal, Exponential, Poisson, Binomial, closed-form MLEs are used when available."),
        tags$li("Gamma MLEs are computed via numerical optimization."),
        tags$li("Contours mark the log-likelihood surface; stars show the MLE; circles show the true parameter.")
      )
    )
  )
)

# ---------- SERVER ----------
server <- function(input, output, session) {
  
  dat <- eventReactive(input$go, {
    set.seed(input$seed)
    n <- input$n
    dist <- input$dist
    x <- switch(dist,
                "Normal"     = rnorm(n, mean = input$norm_mu, sd = input$norm_sigma),
                "Exponential"= rexp(n, rate = input$exp_lambda),
                "Poisson"    = rpois(n, lambda = input$pois_lambda),
                "Binomial"   = rbinom(n, size = input$binom_m, prob = input$binom_p),
                "Gamma"      = rgamma(n, shape = input$gamma_k, scale = input$gamma_theta)
    )
    x
  })
  
  est <- reactive({
    x <- dat()
    req(length(x) > 0)
    dist <- input$dist
    
    if (dist %in% c("Normal", "Exponential", "Poisson", "Binomial")) {
      if (dist == "Binomial") {
        out <- mle_closed(dist, x, m_trials = input$binom_m)
      } else {
        out <- mle_closed(dist, x)
      }
      out$method <- "closed-form"
      return(list(est = out$est, se = out$se, ll = out$ll, method = out$method))
    }
    
    if (dist == "Gamma") {
      # Initials from method of moments
      m <- mean(x); v <- var(x)
      k0 <- max((m^2)/v, 0.1)
      th0 <- max(v/m, 0.1)
      res <- optim(par = c(log(k0), log(th0)),
                   fn = function(p) -ll_gamma(p, x),
                   method = "BFGS", control = list(reltol = 1e-10))
      k_hat <- exp(res$par[1]); th_hat <- exp(res$par[2])
      ll_val <- ll_gamma(c(log(k_hat), log(th_hat)), x)
      # Approx SEs via observed Fisher (inverse Hessian)
      H <- tryCatch(numDeriv::hessian(function(p) -ll_gamma(p, x), c(log(k_hat), log(th_hat))),
                    error = function(e) NULL)
      se <- c(k = NA_real_, theta = NA_real_)
      if (!is.null(H)) {
        V <- tryCatch(solve(H), error = function(e) NULL)
        if (!is.null(V)) {
          # delta method for exp transform: Var(exp(Z)) ~ (exp(Z))^2 Var(Z)
          se <- c(k = k_hat * sqrt(V[1,1]), theta = th_hat * sqrt(V[2,2]))
        }
      }
      return(list(est = c(k = k_hat, theta = th_hat), se = se, ll = ll_val, method = "optim (BFGS)"))
    }
  })
  
  output$est_table <- renderTable({
    req(est())
    dist <- input$dist
    true <- switch(dist,
                   "Normal"     = c(mu = input$norm_mu, sigma = input$norm_sigma),
                   "Exponential"= c(lambda = input$exp_lambda),
                   "Poisson"    = c(lambda = input$pois_lambda),
                   "Binomial"   = c(p = input$binom_p),
                   "Gamma"      = c(k = input$gamma_k, theta = input$gamma_theta)
    )
    e <- est()
    estv <- e$est
    sev  <- e$se
    tab <- data.frame(
      Parameter = names(estv),
      True = as.numeric(true[names(estv)]),
      MLE  = as.numeric(estv),
      SE   = ifelse(is.na(sev), "", sprintf("%.4f", sev)),
      check.names = FALSE
    )
    tab
  })
  
  output$ll_text <- renderText({
    req(est())
    x <- dat(); dist <- input$dist; e <- est()
    
    ll_true <- switch(dist,
                      "Normal"     = ll_norm(c(input$norm_mu, log(input$norm_sigma)), x),
                      "Exponential"= ll_exp(c(log(input$exp_lambda)), x),
                      "Poisson"    = ll_pois(c(log(input$pois_lambda)), x),
                      "Binomial"   = ll_binom(c(qlogis(input$binom_p)), x, input$binom_m),
                      "Gamma"      = ll_gamma(c(log(input$gamma_k), log(input$gamma_theta)), x)
    )
    paste0("Log-likelihood at MLE: ", sprintf("%.3f", e$ll),
           "\nLog-likelihood at True params: ", sprintf("%.3f", ll_true),
           "\nMethod: ", e$method)
  })
  
  output$fitPlot <- renderPlot({
    req(est())
    x <- dat()
    dist <- input$dist
    e <- est()
    
    par(mar = c(4.2, 4.2, 2.5, 1))
    if (dist %in% c("Poisson", "Binomial")) {
      # Discrete: barplot of counts with pmf overlays
      tab <- table(x)
      k <- as.integer(names(tab))
      barplot(tab / length(x), space = 0.2, col = "gray80",
              ylab = "Rel. Frequency", xlab = "Value",
              main = paste0(dist, ": data with pmf (true dashed, MLE solid)"))
      xs <- min(k):max(k)
      
      if (dist == "Poisson") {
        lines(xs - min(xs) + 0.7, dpois(xs, lambda = e$est["lambda"]),
              type = "h", lwd = 3, col = "black")
        lines(xs - min(xs) + 0.7, dpois(xs, lambda = input$pois_lambda),
              type = "h", lwd = 3, col = "red", lty = 2)
      } else { # Binomial
        m <- input$binom_m
        lines(xs - min(xs) + 0.7, dbinom(xs, size = m, prob = e$est["p"]),
              type = "h", lwd = 3, col = "black")
        lines(xs - min(xs) + 0.7, dbinom(xs, size = m, prob = input$binom_p),
              type = "h", lwd = 3, col = "red", lty = 2)
      }
      legend("topright", c("MLE", "True"), lwd = c(3,3), lty = c(1,2),
             col = c("black", "red"), bty = "n")
    } else {
      # Continuous: histogram + densities
      hist(x, freq = FALSE, col = "gray90", border = "white",
           main = paste0(dist, ": data with density (true dashed, MLE solid)"),
           xlab = "x")
      xs <- seq(min(x), max(x), length.out = 400)
      if (dist == "Normal") {
        lines(xs, dnorm(xs, mean = e$est["mu"], sd = e$est["sigma"]), lwd = 3)
        lines(xs, dnorm(xs, mean = input$norm_mu, sd = input$norm_sigma), lwd = 3, col = "red", lty = 2)
      } else if (dist == "Exponential") {
        lines(xs, dexp(xs, rate = e$est["lambda"]), lwd = 3)
        lines(xs, dexp(xs, rate = input$exp_lambda), lwd = 3, col = "red", lty = 2)
      } else if (dist == "Gamma") {
        lines(xs, dgamma(xs, shape = e$est["k"], scale = e$est["theta"]), lwd = 3)
        lines(xs, dgamma(xs, shape = input$gamma_k, scale = input$gamma_theta), lwd = 3, col = "red", lty = 2)
      }
      legend("topright", c("MLE", "True"), lwd = c(3,3), lty = c(1,2),
             col = c("black", "red"), bty = "n")
    }
  })
  
  output$llContour <- renderPlot({
    req(est())
    x <- dat(); dist <- input$dist
    if (!(dist %in% c("Normal", "Gamma"))) return(invisible())
    
    e <- est()
    par(mar = c(4.2, 4.2, 2.5, 1))
    
    if (dist == "Normal") {
      mu_hat <- e$est["mu"]; sig_hat <- e$est["sigma"]
      mu_true <- input$norm_mu; sig_true <- input$norm_sigma
      mu_grid <- seq(mu_hat - 3*sig_hat, mu_hat + 3*sig_hat, length.out = 80)
      sig_grid <- exp(seq(log(sig_hat) - 1.5, log(sig_hat) + 1.5, length.out = 80))
      LL <- outer(mu_grid, sig_grid,
                  Vectorize(function(m,s) ll_norm(c(m, log(s)), x)))
      contour(mu_grid, sig_grid, LL, nlevels = 12,
              xlab = expression(mu), ylab = expression(sigma),
              main = "Normal: log-likelihood contour")
      points(mu_hat, sig_hat, pch = 8, cex = 1.2) # MLE
      points(mu_true, sig_true, pch = 1, col = "red", lwd = 2)
      legend("topright", c("MLE", "True"), pch = c(8,1), col = c("black","red"), bty = "n")
    } else if (dist == "Gamma") {
      k_hat <- e$est["k"]; th_hat <- e$est["theta"]
      k_true <- input$gamma_k; th_true <- input$gamma_theta
      k_grid <- exp(seq(log(k_hat) - 1.2, log(k_hat) + 1.2, length.out = 80))
      th_grid <- exp(seq(log(th_hat) - 1.2, log(th_hat) + 1.2, length.out = 80))
      LL <- outer(k_grid, th_grid,
                  Vectorize(function(k,t) ll_gamma(c(log(k), log(t)), x)))
      contour(k_grid, th_grid, LL, nlevels = 12,
              xlab = expression(k), ylab = expression(theta),
              main = "Gamma: log-likelihood contour")
      points(k_hat, th_hat, pch = 8, cex = 1.2) # MLE
      points(k_true, th_true, pch = 1, col = "red", lwd = 2)
      legend("topright", c("MLE", "True"), pch = c(8,1), col = c("black","red"), bty = "n")
    }
  })
}

shinyApp(ui, server)