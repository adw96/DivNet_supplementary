## @knitr metrics

mse_loss_shannon <- new_metric("mse_loss_shannon", "MSE for estimating Shannon",
                               metric = function(model, out) {
                                 if (is.na(out$shannon) %>% any) {
                                   print(paste("out$shannon:", out$shannon))
                                   stop("oh no! NA for shannon!")
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

true_shannon_list <- new_metric("true_shannon", "true_shannon",
                           metric = function(model, out) {
                             list("shannon" = model$shannon)
                           })


