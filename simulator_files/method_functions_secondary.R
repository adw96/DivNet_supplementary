## @knitr methods

plug_add_half <- new_method("plug-in-add-half", "Imputed Multinomial",
                            method = function(model, draw) {
                              draw_imputed <- draw
                              draw_imputed[which(draw_imputed == 0, arr.ind = TRUE)] <- 0.5
                              draw_imputed <- draw_imputed / apply(draw_imputed, 1, sum)
                              list(shannon = apply(draw_imputed, 1, breakaway::true_shannon),
                                   simpson = apply(draw_imputed, 1, breakaway::true_simpson),
                                   bray_curtis = DivNet::bray_curtis_true(draw_imputed), 
                                   euclidean = DivNet::euclidean_true(draw_imputed))
                            })

## Old

# plug_pooled <- new_method("plug-in-pooled", "Pooled Multinomial",
#                           method = function(model, draw) {
#                             
#                             n <- model$n
#                             if (model$p == 2) {
#                               
#                               shannon1 = vector("numeric", n)
#                               simpson1 = vector("numeric", n)
#                               bray_curtis1 = matrix(0, nrow=n, ncol=n) 
#                               euclidean1 = matrix(0, nrow=n, ncol=n)
#                               
#                               settings <- model$my_x[,2] %>% unique
#                               
#                               for (j in 1:length(settings)) {
#                                 these <- which(model$my_x[,2] == settings[j])
#                                 collapsed <- apply(draw[these, ], 2, sum)
#                                 shannon1[these] <- shannon_true(collapsed)
#                                 simpson1[these] <- simpson_true(collapsed)
#                               }
#                               
#                               if (length(settings) == 2) {
#                                 collapsed1 <- apply(draw[which(model$my_x[,2] == settings[1]), ], 2, sum)
#                                 collapsed2 <- apply(draw[which(model$my_x[,2] == settings[2]), ], 2, sum)
#                                 bc <- bray_curtis2_row(collapsed1, collapsed2)
#                                 euc <- euc_fast(collapsed1, collapsed2)
#                                 for (i in 1:n) {
#                                   j <- which(model$my_x[,2] != model$my_x[i,2])
#                                   bray_curtis1[i,j] <- bc
#                                   euclidean1[i,j] <- euc
#                                 }
#                               }
#                               
#                               list(shannon = shannon1,
#                                    simpson = simpson1,
#                                    bray_curtis = bray_curtis1, 
#                                    euclidean = euclidean1)
#                             } else {
#                               list(shannon = NA,
#                                    simpson = NA,
#                                    bray_curtis = NA, 
#                                    euclidean = NA)
#                             }
#                           })
# 
# miller_maddow <- new_method("miller-maddow", "MLE with Miller Maddow correction",
#                             method = function(model, draw) {
#                               number_species <- apply(draw, 1, function(x) sum(x>0))
#                               ms <- apply(draw, 1, sum)
#                               correction <- (number_species - 1)/(2*ms)
#                               list(shannon = apply(draw, 1, shannon_true) + correction,
#                                    simpson = NA,
#                                    bray_curtis = NA, 
#                                    euclidean = NA)
#                             })
# 
# zhang2012 <- new_method("zhang2012", "Zhang 2012",
#                         method = function(model, draw) {
#                           
#                           zhang <- function(ps, n) {
#                             inner_pieces <- function(v) {
#                               js <- 0:(v-1)
#                               first_piece <- n^(v+1) * factorial(n-v-1) / factorial(n)
#                               second_piece <- sum(ps*prod(1-ps-js/n))
#                               first_piece * second_piece/v
#                             }
#                             sapply(1:(n-1), inner_pieces) %>% sum
#                           }
#                           
#                           ns <- apply(draw, 1, sum)
#                           shannon_hat <- apply(draw, 1, breakaway::make_frequency_count_table) %>% 
#                             lapply(breakaway::to_proportions, type = "frequency count") %>% 
#                             mapply(FUN = zhang, ., n = ns) %>% unlist 
#                           
#                           list(shannon = shannon_hat,
#                                simpson = NA,
#                                bray_curtis = NA, 
#                                euclidean = NA)
#                         })