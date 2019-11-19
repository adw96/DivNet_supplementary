## @knitr methods

divnet_method <- new_method("divnet", "DivNet",
                            method = function(model, draw) {
                              output <- DivNet::divnet(W = draw, 
                                                       X = model$my_x,
                                                       variance="none", 
                                                       base=1)
                              list(shannon = output$shannon,
                                   simpson = output$simpson,
                                   bray_curtis = output[["bray-curtis"]], 
                                   euclidean = output$euclidean
                              )
                            })

plug <- new_method("plug-in", "Multinomial",
                   method = function(model, draw) {
                     
                     rel_abunds <- draw / (draw %>% apply(1, sum))
                     list(shannon = apply(draw, 1, breakaway::sample_shannon) %>% 
                            lapply(function(x) {x$estimate}) %>% unlist,
                          simpson = apply(draw, 1, breakaway::sample_simpson) %>% 
                                            lapply(function(x) {x$estimate}) %>% unlist,
                          bray_curtis = DivNet::bray_curtis_true(rel_abunds), 
                          euclidean = DivNet::euclidean_true(rel_abunds))
                   })


chao_shen <- new_method("chao-shen", "Chao-Shen",
                        method = function(model, draw) {
                          f_tables <- apply(draw, 1, breakaway::make_frequency_count_table)
                          f_tables2 <- f_tables %>% lapply(breakaway:::check_format)
                          ests <- f_tables2 %>% lapply(breakaway:::chao_shen_estimate) %>% unlist
                          
                          list(shannon = ests,
                               simpson = NA,
                               bray_curtis = NA, 
                               euclidean = NA)
                        })


inext <- new_method("inext", "iNEXT",
                    method = function(model, draw) {
                      reformat_draw <- draw %>% t
                      rownames(reformat_draw) <- paste("Species", 1:(dim(reformat_draw)[1]), sep="")
                      colnames(reformat_draw) <- paste("Site", 1:(dim(reformat_draw)[2]), sep="")
                      my_table <- iNEXT::iNEXT(reformat_draw)$AsyEst
                      list(shannon = my_table[my_table[, "Diversity"] ==  "Shannon diversity", "Estimator"] %>% log, 
                           simpson = 1/my_table[my_table[, "Diversity"] == "Simpson diversity", "Estimator"],
                           bray_curtis = NA, 
                           euclidean = NA)
                    })

arbel <- new_method("arbel", "Arbel et. al", 
                    method = function(model, draw) {
                      
                      p <- dim(model$my_x)[2]
                      
                      if (p == 2) {
                        
                        get_shannon <- function(my_matrix) {
                          my_matrix %>% apply(1, DivNet::shannon_true)
                        }
                        
                        get_simpson <- function(my_matrix) {
                          my_matrix %>% apply(1, DivNet::simpson_true)
                        }
                        
                        output <- DepGEM::gibbs(n.iter=500, Y = draw, X = model$my_x[,2])
                        
                        zz <- output$p_store 
                        
                        list(shannon = apply(zz, 3, get_shannon) %>% apply(1, mean),
                             simpson = apply(zz, 3, get_simpson) %>% apply(1, mean),
                             bray_curtis = DivNet::bray_curtis_true(apply(zz, c(1,2), mean)), 
                             euclidean = DivNet::euclidean_true(apply(zz, c(1,2), mean)))
                      } else {
                        list(shannon = NA,
                             simpson = NA,
                             bray_curtis = NA, 
                             euclidean = NA)
                      }
                    })