## @knitr models

simulate_data_sigma <- function(q, beta_sd) {
  
  # find the q/2 most abundant taxa and the q/2 randomly sampled ones
  random_sample <- sample(1:length(names(my_w)), floor(q/2))
  most_abundant <- apply(my_w[, -random_sample], 2, sum) %>% sort(decreasing = TRUE) %>% head(ceiling(q/2)) %>% names
  chosen <- c(names(my_w)[random_sample], most_abundant)
  
  # let's get all of the altered and glassy basalts and subset to the taxa we chose
  W <- my_w[-11, chosen]
  X <- my_x[-11, -2]
  ms <- apply(my_w, 1, sum)[-11]
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  beta <- matrix(rnorm(p * (q-1), mean=0, sd = beta_sd), nrow=p, ncol=q-1) # p x q
  mu <- X %*% beta
  
  ## set up covariance
  Y.p <- to_log_ratios(W, base = 1, p = 0.5)
  # edit following first revisions: now have dense matrix
  Sigma <- cov(Y.p)
  
  new_model(name = "aitchison-model",
            label = sprintf("Lee data (q = %s, n = %s, p = 2)", 
                            q,
                            n),
            params = list(W = W,
                          X = X,
                          q = q,
                          n = n,
                          p = p,
                          mu = mu, 
                          Sigma = Sigma,
                          ms = ms
                          ),
            simulate = function(mu, Sigma, nsim) {
              # this function must return a list of length nsim
              replicate(nsim,
                        make_w(mu, Sigma, mm = ms, base=1), simplify=F)
            }
  )
}


simulate_data_sigma_rich_structure <- function(q, beta_sd) {
  
  # find the q most abundant taxa 
  most_abundant <- apply(my_w, 2, sum) %>% sort(decreasing = TRUE) %>% head(ceiling(q)) %>% names

  # let's get all of the altered and glassy basalts and subset to the taxa we chose
  W <- my_w[-11, most_abundant]
  X <- my_x[-11, -2]
  ms <- apply(my_w, 1, sum)[-11]
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  Y.p <- to_log_ratios(W, base = 1, p = 0.5)
  Sigma <- cov(Y.p)
  beta <- DivNet:::OLS(X, Y.p)
  mu <- X %*% beta
  
  new_model(name = "aitchison-model",
            label = sprintf("Lee data (q = %s, n = %s, p = 2)", 
                            q,
                            n),
            params = list(W = W,
                          X = X,
                          q = q,
                          n = n,
                          p = p,
                          mu = mu, 
                          Sigma = Sigma,
                          ms = ms
            ),
            simulate = function(mu, Sigma, nsim) {
              # this function must return a list of length nsim
              replicate(nsim,
                        make_w(mu, Sigma, mm = ms, base=1), simplify=F)
            }
  )
}
