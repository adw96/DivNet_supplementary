## @knitr methods

divnet_method_lv <- new_method("divnet-lv", "DivNet (LV)",
                               method = function(model, draw) {
                                 out <- DivNet::divnet(W = draw, 
                                                          X = matrix(1, ncol = 1, nrow = model$ll + 1), 
                                                          variance="none", 
                                                          base=1)
                                 
                                 shann_ests <- out$shannon %>% 
                                   lapply(function(x) {x$estimate}) %>% 
                                   unlist %>% mean
                                 
                                 list(shannon = shann_ests)
                               })

plug_lv <- new_method("plug-in-lv", "Multinomial (LV)",
                      method = function(model, draw) {
                        
                        rel_abunds <- draw / (draw %>% apply(1, sum))
                        
                        shann_ests <- apply(draw, 1, breakaway::sample_shannon) %>% 
                          lapply(function(x) {x$estimate}) %>% 
                          unlist %>% mean
                        
                        list(shannon = shann_ests)
                      })

chao_shen_lv <- new_method("chao-shen-lv", "Chao-Shen (LV)",
                           method = function(model, draw) {
                             f_tables <- apply(draw, 1, breakaway::make_frequency_count_table)
                             f_tables2 <- f_tables %>% lapply(breakaway:::check_format)
                             ests <- f_tables2 %>% lapply(breakaway:::chao_shen_estimate)
                             
                             shann_ests <- ests %>% 
                               #lapply(function(x) {x$estimate}) %>% 
                               unlist %>% mean
                             
                             list(shannon = shann_ests)
                           })

inext_lv <- new_method("inext-lv", "iNEXT (LV)",
                    method = function(model, draw) {
                      reformat_draw <- draw %>% t
                      rownames(reformat_draw) <- paste("Species", 1:(dim(reformat_draw)[1]), sep="")
                      colnames(reformat_draw) <- paste("Site", 1:(dim(reformat_draw)[2]), sep="")
                      my_table <- iNEXT::iNEXT(reformat_draw)$AsyEst
                      list(shannon = my_table[my_table[, "Diversity"] ==  "Shannon diversity", "Estimator"] %>% log %>% mean)
                    })

arbel_lv <- new_method("arbel-lv", "Arbel et. al (LV)", 
                       method = function(model, draw) {
                         
                         get_shannon <- function(my_matrix) {
                           my_matrix %>% apply(1, breakaway::true_shannon)
                         }
                         
                         output <- DepGEM::gibbs(n.iter=500, 
                                                 Y = draw, 
                                                 X = 1)
                         
                         zz <- output$p_store 
                         
                         stop("no")
                         list(shannon = apply(zz, 3, get_shannon) %>% apply(1, mean))
                       })
