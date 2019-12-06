## @knitr models
nonlinear_function1 <- function(x, q, gamma1) {
  1-gamma1*x*(x-10)
}

make_aitchison_nonlinear_function1 <- function(n, q, 
                                               gamma1, 
                                               sigma_max, sigma_min, 
                                               mm = 1e5) {
  
  # make covariance matrix
  A <- matrix(runif((q-1)^2, min=-1, max=1), ncol=q - 1) 
  D <- diag(seq(from = sigma_max, to=sigma_min, length.out=q - 1))
  Sigma <- t(A) %*% D %*% A
  
  # make mu 
  xx <- seq(0, 10, length.out=n)
  mu <- matrix(0, nrow = n, ncol = q-1)
  mu[, 1] <- nonlinear_function1(x = xx, q = q, gamma1 = gamma1)
  mu_props <- DivNet::to_composition_matrix(Y=mu, base=1)
  
  my_x <- cbind(rep(1, n), xx - mean(xx), (xx - mean(xx))^2)
  
  new_model(name = "nonlinear-aitchison",
            label = sprintf("Nonlinear (n = %s, q = %s, gamma1 = %s)", n, q, gamma1),
            params = list(n = n, 
                          q = q, 
                          mu = mu, 
                          Sigma = Sigma,
                          xx = xx,
                          my_x = my_x, 
                          shannon = apply(mu_props, 1, breakaway::true_shannon),
                          mm = mm),
            simulate = function(mu, Sigma, mm, q, nsim) {
              # this function must return a list of length nsim
              replicate(nsim,
                        DivNet::make_w(mu, Sigma, mm=mm, base=1), simplify=F)
            })
}

nonlinear_function2 <- function(x, q, gamma1) {
  3 + 0.25*exp(gamma1*x)
}

make_aitchison_nonlinear_function2 <- function(n, q, 
                                               gamma1, 
                                               sigma_max, sigma_min, 
                                               mm = 1e5) {
  
  # make covariance matrix
  A <- matrix(runif((q-1)^2, min=-1, max=1), ncol=q - 1) 
  D <- diag(seq(from = sigma_max, to=sigma_min, length.out=q - 1))
  Sigma <- t(A) %*% D %*% A
  
  # make mu 
  xx <- seq(0, 10, length.out=n)
  mu <- matrix(0, nrow = n, ncol = q-1)
  mu[, 1] <- nonlinear_function2(x = xx, q = q, gamma1 = gamma1)
  mu_props <- DivNet::to_composition_matrix(Y=mu, base=1)
  
  my_x <- cbind(rep(1, n), xx - mean(xx), (xx - mean(xx))^2)
  
  new_model(name = "nonlinear-aitchison-2",
            label = sprintf("Nonlinear Model 2(n = %s, q = %s, gamma1 = %s)", n, q, gamma1),
            params = list(n = n, 
                          q = q, 
                          mu = mu, 
                          Sigma = Sigma,
                          xx = xx,
                          my_x = my_x, 
                          shannon = apply(mu_props, 1, breakaway::true_shannon),
                          mm = mm),
            simulate = function(mu, Sigma, mm, q, nsim) {
              # this function must return a list of length nsim
              replicate(nsim,
                        DivNet::make_w(mu, Sigma, mm=mm, base=1), simplify=F)
            })
}
