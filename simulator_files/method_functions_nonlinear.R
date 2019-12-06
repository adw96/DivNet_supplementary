## @knitr methods

divnet_method_nonlinear <- new_method("divnet-nonlin", "DivNet (Nonlinear)",
                                      method = function(model, draw) {
                                        output <- DivNet::divnet(W = draw, 
                                                                 X = model$my_x,
                                                                 variance="none", 
                                                                 base=2)
                                        list(shannon = output$shannon)
                                      })
divnet_method_linear <- new_method("divnet-lin", "DivNet (Linear)",
                                   method = function(model, draw) {
                                     output <- DivNet::divnet(W = draw, 
                                                              X = model$my_x[,1:2],
                                                              variance="none", 
                                                              base=2)
                                     list(shannon = output$shannon)
                                   })

plug <- new_method("plug-in", "Multinomial",
                   method = function(model, draw) {
                     rel_abunds <- draw / (draw %>% apply(1, sum))
                     list(shannon = apply(draw, 1, breakaway::sample_shannon) %>% 
                            lapply(function(x) {x$estimate}) %>% unlist)
                   })


chao_shen <- new_method("chao-shen", "Chao-Shen",
                        method = function(model, draw) {
                          f_tables <- apply(draw, 1, breakaway::make_frequency_count_table)
                          f_tables2 <- f_tables %>% lapply(breakaway:::check_format)
                          ests <- f_tables2 %>% lapply(breakaway:::chao_shen_estimate) %>% unlist
                          
                          list(shannon = ests)
                        })


inext <- new_method("inext", "iNEXT",
                    method = function(model, draw) {
                      reformat_draw <- draw %>% t
                      rownames(reformat_draw) <- paste("Species", 1:(dim(reformat_draw)[1]), sep="")
                      colnames(reformat_draw) <- paste("Site", 1:(dim(reformat_draw)[2]), sep="")
                      my_table <- iNEXT::iNEXT(reformat_draw)$AsyEst
                      list(shannon = my_table[my_table[, "Diversity"] ==  "Shannon diversity", "Estimator"] %>% log)
                    })

arbel_linear <- new_method("arbel-linear", "Arbel et. al (Linear)", 
                           method = function(model, draw) {
                             
                             get_shannon <- function(my_matrix) {
                               my_matrix %>% apply(1, breakaway::true_shannon)
                             }
                             
                             output <- DepGEM::gibbs(n.iter=500, Y = draw, 
                                                     X = model$my_x[,2])
                             
                             zz <- output$p_store 
                             
                             list(shannon = apply(zz, 3, get_shannon) %>% apply(1, mean))
                             
                           })

arbel_quad <- new_method("arbel-quad", "Arbel et. al (Quadratic)", 
                           method = function(model, draw) {
                             
                             get_shannon <- function(my_matrix) {
                               my_matrix %>% apply(1, breakaway::true_shannon)
                             }
                             
                             output <- DepGEM::gibbs(n.iter=500, 
                                                     Y = draw, 
                                                     X = model$my_x[,3])
                             
                             zz <- output$p_store 
                             
                             list(shannon = apply(zz, 3, get_shannon) %>% apply(1, mean))
                             
                           })
