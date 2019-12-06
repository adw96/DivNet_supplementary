n = 20 
q = 20 
gamma1 = 0.5 
sigma_max = 5 
sigma_min = 0.01
mm = 1e5

A <- matrix(runif((q-1)^2, min=-1, max=1), ncol=q - 1) 
D <- diag(seq(from = sigma_max, to=sigma_min, length.out = q - 1))
Sigma <- t(A) %*% D %*% A

set.seed(1)
xx <- seq(0, 10, length.out=n)
mu <- matrix(0, nrow = n, ncol = q-1)
mu[, 1] <- nonlinear_function1(x = xx, q = q, gamma1 = gamma1)
mu_props <- DivNet::to_composition_matrix(Y=mu, base=1)
true_shannon <- apply(mu_props, 1, breakaway::true_shannon)
ww <- DivNet::make_w(mu, Sigma, mm=mm, base=1)
ww
output <- DivNet::divnet(W = ww, 
                         X = my_x,
                         variance="none", 
                         base=2)
# 
# mu_props
# mu_props %>% dim
# 
# my_x <- cbind(rep(1, n), xx - mean(xx), (xx - mean(xx))^2)
# 
# true_shannon <- apply(mu_props, 1, breakaway::true_shannon)
# ww <- DivNet::make_w(mu, Sigma, mm=mm, base=q)
# ww
output <- DivNet::divnet(W = ww, 
                         X = my_x,
                         variance="none", 
                         base=2)
output
data.frame(x = xx, 
           fitted_z = output$fitted_z[,1], 
           mu_props = mu_props[, 1]) %>%
  ggplot(aes(x, fitted_z)) +
  geom_point(col = "red") + 
  geom_point(aes(y = mu_props)) 

data.frame(x = xx, 
           fitted_z = output$fitted_z[,2], 
           mu_props = mu_props[, 2]) %>%
  ggplot(aes(x, fitted_z)) +
  geom_point(col = "red") + 
  geom_point(aes(y = mu_props)) 

data.frame(x = xx, 
           true_shannon = true_shannon, 
           estimated_shannon = output$shannon %>% 
             lapply(function(x) {x$estimate}) %>% unlist) %>%
  ggplot(aes(x, true_shannon)) +
  geom_point(col = "red") + 
  geom_point(aes(y = estimated_shannon)) 
