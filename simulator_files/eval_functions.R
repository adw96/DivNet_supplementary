## @knitr metrics

mse_loss_shannon <- new_metric("mse_loss_shannon", "MSE for estimating Shannon",
                               metric = function(model, out) {
                                 if (is.na(out$shannon) %>% any) {
                                   NA
                                 } else if ((model$shannon %>% length) != (out$shannon %>% length)) {
                                   print(paste("model$shannon:", model$shannon))
                                   print(paste("out$shannon:", out$shannon))
                                   stop("Wrong lengths for shannon index")  
                                 } else if ("alpha_estimates" %in% class(out$shannon)) {
                                   
                                   shann_ests <- out$shannon %>% 
                                     lapply(function(x) {x$estimate}) %>% 
                                     unlist
                                   
                                   (model$shannon - shann_ests)^2 %>% mean
                                   
                                 } else {
                                   (model$shannon - out$shannon)^2 %>% mean
                                 }
                               })

estimated_shannon <- new_metric("estimated_shannon", "estimated_shannon",
                                metric = function(model, out) {
                                  if (is.na(out$shannon) %>% any) {
                                    NA
                                  }
                                  else if ((model$shannon %>% length) != (out$shannon %>% length)) {
                                    print(paste("model$shannon:", model$shannon))
                                    print(paste("out$shannon:", out$shannon))
                                    stop("Wrong lengths for shannon index")  
                                  } else {
                                    out$shannon %>% mean
                                  }
                                })

true_shannon <- new_metric("true_shannon", "true_shannon",
                           metric = function(model, out) {
                             if (is.na(out$shannon) %>% any) {
                               NA
                             }
                             else if ((model$shannon %>% length) != (out$shannon %>% length)) {
                               print(paste("model$shannon:", model$shannon))
                               print(paste("out$shannon:", out$shannon))
                               stop("Wrong lengths for shannon index")  
                             } else {
                               model$shannon %>% mean
                             }
                           })



mse_loss_simpson <- new_metric("mse_loss_simpson", "MSE for estimating Simpson",
                               metric = function(model, out) {
                                 if (is.na(out$simpson) %>% any) {
                                   NA
                                 } else if ((model$simpson %>% length) != (out$simpson %>% length)) {
                                   print(paste("model$simpson:", model$simpson))
                                   print(paste("out$simpson:", out$simpson))
                                   stop("Wrong lengths for simpson index")  
                                 } else if ("alpha_estimates" %in% class(out$simpson)) {
                                   
                                   simp_ests <- out$simpson %>% 
                                     lapply(function(x) {x$estimate}) %>% 
                                     unlist
                                   
                                   (model$simpson - simp_ests)^2 %>% mean
                                 } else {
                                   (model$simpson - out$simpson)^2 %>% mean
                                 }
                               })

mse_loss_bray_curtis <- new_metric("mse_loss_bray_curtis", "MSE for estimating Bray Curtis",
                                   metric = function(model, out) {
                                     (model$bray_curtis - out$bray_curtis)^2 %>% mean
                                   })
bray_curtis_true <- new_metric("bray_curtis_true", "bray_curtis_true",
                               metric = function(model, out) {
                                 (model$bray_curtis) %>% mean
                               })
bray_curtis_est <- new_metric("bray_curtis_est", "bray_curtis_est",
                              metric = function(model, out) {
                                (out$bray_curtis) %>% mean
                              })

mse_loss_euclidean <- new_metric("mse_loss_euclidean", "MSE for estimating Euclidean",
                                 metric = function(model, out) {
                                   (model$euclidean - out$euclidean)^2 %>% mean
                                 })

mse_loss_bray_curtis_subset <- new_metric("mse_loss_bray_curtis_subset", "MSE for estimating Bray Curtis",
                                          metric = function(model, out) {
                                            my_fitted <- model$bray_curtis
                                            my_true <- out$bray_curtis
                                            my_fitted[lower.tri(my_fitted)] <- 0
                                            my_true[lower.tri(my_true)] <- 0
                                            (my_fitted - my_true)^2 %>% mean
                                          })

mse_loss_euclidean_subset <- new_metric("mse_loss_euclidean_subset", "MSE for estimating Euclidean",
                                        metric = function(model, out) {
                                          my_fitted <- model$euclidean
                                          my_true <- out$euclidean
                                          my_fitted[lower.tri(my_fitted)] <- 0
                                          my_true[lower.tri(my_true)] <- 0
                                          (my_fitted - my_true)^2 %>% mean
                                        })




