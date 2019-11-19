## @knitr models

make_aitchison <- function(n, p, q, 
                           sigma_max, sigma_min, beta_sd,
                           mm = 1e5) {
  
  # make covariance matrix
  A <- matrix(runif(q^2, min=-1, max=1), ncol=q) 
  D <- diag(seq(from = sigma_max, to=sigma_min, length.out=q))
  Sigma <- t(A) %*% D %*% A
  
  # make mu
  my_x <- cbind(rep(1, n))
  beta <- matrix(rnorm(p * q, mean=0, sd = beta_sd), nrow=p, ncol=q) # p x q
  if (p == 1) {
    my_x <- cbind(rep(1, n)) # n x 1
  } else if (p  == 2) {
    stopifnot(n %% 2 == 0)
    my_x <- cbind(rep(1, n), rep(c(0,1), each = n/2)) # n x 2
  } else if (p == 3) {
    stopifnot(n %% 2 == 0)
    my_x <- cbind(rep(1, n), rep(c(0,1), each = n/2), rep(c(0,1), n/2)) # n x 3
  }
  mu <- my_x %*% beta 
  
  mu_props <-  DivNet::to_composition_matrix(Y=mu, base=1)
  bc_true <- bray_curtis_true(mu_props)
  euclidean_true <- euclidean_true(mu_props)
  
  new_model(name = "aitchison-model-alpha",
            label = sprintf("Aitchison (n = %s, q = %s, p = %s, sigma_min = %s, sigma_max = %s)", n, q, p, sigma_min, sigma_max),
            params = list(n = n, 
                          mu = mu, 
                          Sigma = Sigma,
                          sigma_max = sigma_max, 
                          sigma_min = sigma_min, 
                          p = p, 
                          q = q, 
                          beta = beta,
                          my_x = my_x,
                          shannon = apply(mu_props, 1, breakaway::true_shannon),
                          simpson = apply(mu_props, 1, breakaway::true_simpson),
                          bray_curtis = bc_true,
                          euclidean = euclidean_true,
                          mm = mm),
            simulate = function(mu, Sigma, mm, nsim) {
              # this function must return a list of length nsim
              replicate(nsim,
                        DivNet::make_w(mu, Sigma, mm=mm, base=1), simplify=F)
              
            })
}


make_fisher_mehta <- function(n, q) {
  
  q = 10
  n = 20
  eta <- rnorm(q*n, mean=0, sd=0.1) %>% exp %>% matrix(., nrow=n, ncol=q)
  x_steady <- rnorm(q, mean=1, sd = 0.1) %>% exp
  x <- matrix(NA, nrow=n, ncol = q)
  x[1, ] <- 1
  delta <- 0.1
  cc <- matrix(0, ncol = q, nrow = q)
  diag(cc) <- runif(n=q, min = -1.9, max = -0.1) / x_steady
  possible_indices <- which(cc == 0, arr.ind=TRUE)
  selected_offdiagonals <- sample(1:nrow(possible_indices), size = 10, replace=FALSE)
  off_diagonal_values <- runif(n=10, min = -1.9, max = -0.1) / x_steady ## Not specified in paper!!!
  
  possible_indices[selected_diagonals, ]
  
  x_steady <- rnorm(q, mean=1, sd = 0.1) %>% exp
  
  
  # make covariance matrix
  A <- matrix(runif(q^2, min=-1, max=1), ncol=q) 
  D <- diag(seq(from = sigma_max, to=sigma_min, length.out=q))
  Sigma <- t(A) %*% D %*% A
  
  # make mu
  my_x <- cbind(rep(1, n))
  beta <- matrix(rnorm(p * q, mean=0, sd = beta_sd), nrow=p, ncol=q) # p x q
  if (p == 1) {
    my_x <- cbind(rep(1, n)) # n x 1
  } else if (p  == 2) {
    stopifnot(n %% 2 == 0)
    my_x <- cbind(rep(1, n), rep(c(0,1), each = n/2)) # n x 2
  } else if (p == 3) {
    stopifnot(n %% 2 == 0)
    my_x <- cbind(rep(1, n), rep(c(0,1), each = n/2), rep(c(0,1), n/2)) # n x 3
  }
  mu <- my_x %*% beta 
  
  mu_props <-  DivNet::to_composition_matrix(Y=mu, base=1)
  bc_true <- bray_curtis_true(mu_props)
  euclidean_true <- euclidean_true(mu_props)
  
  new_model(name = "aitchison-model-alpha",
            label = sprintf("Aitchison (n = %s, q = %s, p = %s, sigma_min = %s, sigma_max = %s)", n, q, p, sigma_min, sigma_max),
            params = list(n = n, 
                          mu = mu, 
                          Sigma = Sigma,
                          sigma_max = sigma_max, 
                          sigma_min = sigma_min, 
                          p = p, 
                          q = q, 
                          beta = beta,
                          my_x = my_x,
                          shannon = apply(mu_props, 1, breakaway::true_shannon),
                          simpson = apply(mu_props, 1, breakaway::true_simpson),
                          bray_curtis = bc_true,
                          euclidean = euclidean_true,
                          mm = mm),
            simulate = function(mu, Sigma, mm, nsim) {
              # this function must return a list of length nsim
              replicate(nsim,
                        DivNet::make_w(mu, Sigma, mm=mm), simplify=F)
              
            })
}