##############################################################################################################
##### play around with different covariance measurements and different methods of generating variance #####
##############################################################################################################

make_parametric_variance_estimate <- function(nn) {
  if (nn == 0) {
    cov_method <- "default"
  } else if (nn == 1) {
    cov_method <- "diagonal"
  } else if (nn == 2) {
    cov_method <- "stars"
  } else {
    stop("Unsure of the variance estimation procedure")
  }
  new_method(name=sprintf("parametric_%s",  cov_method), 
             label = sprintf("Parametric, %s",  cov_method),
             settings=list(cov_method = cov_method),
             method = function(model, draw, cov_method) {
               
               dv = divnet(draw, model$X, tuning="fast", 
                           network=cov_method, 
                           B = 3, 
                           variance="parametric")
               
               shannon_variance_estimated <- dv$`shannon-variance`
               bc_variance_estimated <- dv$`bray-curtis-variance`
               
               #### true variance
               simulated_data <- replicate(3,
                                           make_w(model$mu, model$Sigma, model$ms), simplify=F)
               simulated_data <- replicate(3,
                                           make_w(mu, Sigma, ms), simplify=F)
               
               estimated_diversity <- lapply(simulated_data, 
                                             FUN=divnet, 
                                             model$X, tuning="fast", 
                                             network=cov_method, variance = "none") 
               
               shannon_variance_true <- estimated_diversity %>% 
                 lapply(function(x) x$fitted_z) %>%
                 lapply(function(x) apply(x, 1, shannon)) %>%
                 simplify2array %>% 
                 apply(1, var)
               
               bc_variance_true <- estimated_diversity %>% 
                 lapply(function(x) x$fitted_z) %>%
                 lapply(function(x) bray_curtis_true(x)) %>% 
                 simplify2array %>% apply(1:2, var)
               
               list("shannon_variance_estimated" = shannon_variance_estimated,
                    "shannon_variance_true" = shannon_variance_true,
                    "bc_variance_estimated" = bc_variance_estimated,
                    "bc_variance_true" = bc_variance_true)
             }
  )
}

make_nonparametric_variance_estimate <- function(nn) {
  if (nn == 0) {
    cov_method <- "default"
  } else if (nn == 1) {
    cov_method <- "diagonal"
  } else if (nn == 2) {
    cov_method <- "stars"
  } else {
    stop("Unsure of the variance estimation procedure")
  }
  new_method(name=sprintf("nonparametric_%s",  cov_method), 
             label = sprintf("Nonparametric, %s",  cov_method),
             settings=list(cov_method = cov_method),
             
             method = function(model, draw, cov_method) {
               
               dv = divnet(draw, model$X, tuning="fast", 
                           network=cov_method, 
                           B = 3, 
                           variance="nonparametric")
               
               shannon_variance_estimated <- dv$`shannon-variance`
               bc_variance_estimated <- dv$`bray-curtis-variance`
               
               #### true variance
               simulated_data <- replicate(3,
                                           make_w(model$mu, model$Sigma, model$ms), simplify=F)
               
               estimated_diversity <- lapply(simulated_data, 
                                             FUN=divnet, 
                                             model$X, tuning="fast", 
                                             network=cov_method, variance = "none") 
               
               shannon_variance_true <- estimated_diversity %>% 
                 lapply(function(x) x$fitted_z) %>%
                 lapply(function(x) apply(x, 1, shannon)) %>%
                 simplify2array %>% 
                 apply(1, var)
               
               bc_variance_true <- estimated_diversity %>% 
                 lapply(function(x) x$fitted_z) %>%
                 lapply(function(x) bray_curtis_true(x)) %>% 
                 simplify2array %>% 
                 apply(1:2, var)
               
               print("bc_variance_true")
               print(bc_variance_true)
               print("bc_variance_estimated")
               print(bc_variance_estimated)

               list("shannon_variance_estimated" = shannon_variance_estimated,
                    "shannon_variance_true" = shannon_variance_true,
                    "bc_variance_estimated" = bc_variance_estimated,
                    "bc_variance_true" = bc_variance_true)
             }
  )
}
