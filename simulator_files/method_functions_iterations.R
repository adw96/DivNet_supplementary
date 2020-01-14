## @knitr methods


divnet_6 <- new_method("divnet6", "DivNet 6",
                       method = function(model, draw) {
                         output <- DivNet::divnet(W = draw, 
                                                  X = model$my_x,
                                                  variance="none", 
                                                  base=1, 
                                                  tuning="fast")
                         list(shannon = output$shannon %>% summary %$% estimate)
                       })


divnet_20 <- new_method("divnet20", "DivNet 20",
                       method = function(model, draw) {
                         output <- DivNet::divnet(W = draw, 
                                                  X = model$my_x,
                                                  variance="none", 
                                                  base=1, 
                                                  tuning=list("EMiter" = 20,
                                                              "EMburn" = 10,
                                                              "MCiter" = 500,
                                                              "MCburn" = 250))
                         list(shannon = output$shannon %>% summary %$% estimate)
                       })


divnet_10 <- new_method("divnet10", "DivNet 10",
                        method = function(model, draw) {
                          output <- DivNet::divnet(W = draw, 
                                                   X = model$my_x,
                                                   variance="none", 
                                                   base=1, 
                                                   tuning="careful")
                          list(shannon = output$shannon %>% summary %$% estimate)
                        })


divnet_50 <- new_method("divnet50", "DivNet 50",
                        method = function(model, draw) {
                          output <- DivNet::divnet(W = draw, 
                                                   X = model$my_x,
                                                   variance="none", 
                                                   base=1, 
                                                   tuning=list("EMiter" = 50,
                                                               "EMburn" = 25,
                                                               "MCiter" = 500,
                                                               "MCburn" = 250))
                          list(shannon = output$shannon %>% summary %$% estimate)
                        })



arbel_100 <- new_method("arbel100", "Arbel et. al 100", 
                        method = function(model, draw) {
                          
                          p <- dim(model$my_x)[2]
                          
                          if (p == 2) {
                            
                            get_shannon <- function(my_matrix) {
                              my_matrix %>% apply(1, breakaway::true_shannon)
                            }
                            
                            get_simpson <- function(my_matrix) {
                              my_matrix %>% apply(1, breakaway::true_simpson)
                            }
                            
                            output <- DepGEM::gibbs(n.iter=100, 
                                                    Y = draw, 
                                                    X = model$my_x[,2])
                            
                            zz <- output$p_store 
                            
                            list(shannon = apply(zz, 3, get_shannon) %>% apply(1, mean),
                                 simpson = apply(zz, 3, get_simpson) %>% apply(1, mean),
                                 bray_curtis = DivNet::bray_curtis_true(apply(zz, c(1,2), mean)), 
                                 euclidean = DivNet::euclidean_true(apply(zz, c(1,2), mean)))
                          } else {
                            stop("p must be 2")
                          }
                        })


arbel_500 <- new_method("arbel500", "Arbel et. al 500", 
                        method = function(model, draw) {
                          
                          p <- dim(model$my_x)[2]
                          
                          if (p == 2) {
                            
                            get_shannon <- function(my_matrix) {
                              my_matrix %>% apply(1, breakaway::true_shannon)
                            }
                            
                            get_simpson <- function(my_matrix) {
                              my_matrix %>% apply(1, breakaway::true_simpson)
                            }
                            
                            output <- DepGEM::gibbs(n.iter=500, 
                                                    Y = draw, 
                                                    X = model$my_x[,2])
                            
                            zz <- output$p_store 
                            
                            list(shannon = apply(zz, 3, get_shannon) %>% apply(1, mean),
                                 simpson = apply(zz, 3, get_simpson) %>% apply(1, mean),
                                 bray_curtis = DivNet::bray_curtis_true(apply(zz, c(1,2), mean)), 
                                 euclidean = DivNet::euclidean_true(apply(zz, c(1,2), mean)))
                          } else {
                            stop("p must be 2")
                          }
                        })

arbel_1000 <- new_method("arbel1000", "Arbel et. al 1000", 
                       method = function(model, draw) {
                         
                         p <- dim(model$my_x)[2]
                         
                         if (p == 2) {
                           
                           get_shannon <- function(my_matrix) {
                             my_matrix %>% apply(1, breakaway::true_shannon)
                           }
                           
                           get_simpson <- function(my_matrix) {
                             my_matrix %>% apply(1, breakaway::true_simpson)
                           }
                           
                           output <- DepGEM::gibbs(n.iter=1000, 
                                                   Y = draw, 
                                                   X = model$my_x[,2])
                           
                           zz <- output$p_store 
                           
                           list(shannon = apply(zz, 3, get_shannon) %>% apply(1, mean),
                                simpson = apply(zz, 3, get_simpson) %>% apply(1, mean),
                                bray_curtis = DivNet::bray_curtis_true(apply(zz, c(1,2), mean)), 
                                euclidean = DivNet::euclidean_true(apply(zz, c(1,2), mean)))
                         } else {
                           stop("p must be 2")
                         }
                       })

arbel_2000 <- new_method("arbel2000", "Arbel et. al 2000", 
                         method = function(model, draw) {
                           
                           p <- dim(model$my_x)[2]
                           
                           if (p == 2) {
                             
                             get_shannon <- function(my_matrix) {
                               my_matrix %>% apply(1, breakaway::true_shannon)
                             }
                             
                             get_simpson <- function(my_matrix) {
                               my_matrix %>% apply(1, breakaway::true_simpson)
                             }
                             
                             output <- DepGEM::gibbs(n.iter=2000, 
                                                     Y = draw, 
                                                     X = model$my_x[,2])
                             
                             zz <- output$p_store 
                             
                             list(shannon = apply(zz, 3, get_shannon) %>% apply(1, mean),
                                  simpson = apply(zz, 3, get_simpson) %>% apply(1, mean),
                                  bray_curtis = DivNet::bray_curtis_true(apply(zz, c(1,2), mean)), 
                                  euclidean = DivNet::euclidean_true(apply(zz, c(1,2), mean)))
                           } else {
                             stop("p must be 2")
                           }
                         })


arbel_250 <- new_method("arbel250", "Arbel et. al 250", 
                         method = function(model, draw) {
                           
                           p <- dim(model$my_x)[2]
                           
                           if (p == 2) {
                             
                             get_shannon <- function(my_matrix) {
                               my_matrix %>% apply(1, breakaway::true_shannon)
                             }
                             
                             get_simpson <- function(my_matrix) {
                               my_matrix %>% apply(1, breakaway::true_simpson)
                             }
                             
                             output <- DepGEM::gibbs(n.iter=250, 
                                                     Y = draw, 
                                                     X = model$my_x[,2])
                             
                             zz <- output$p_store 
                             
                             list(shannon = apply(zz, 3, get_shannon) %>% apply(1, mean),
                                  simpson = apply(zz, 3, get_simpson) %>% apply(1, mean),
                                  bray_curtis = DivNet::bray_curtis_true(apply(zz, c(1,2), mean)), 
                                  euclidean = DivNet::euclidean_true(apply(zz, c(1,2), mean)))
                           } else {
                             stop("p must be 2")
                           }
                         })