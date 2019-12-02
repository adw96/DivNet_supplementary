# This is the main simulator file
directory <- "/Users/adwillis/research/alpha-aitchison/covariance-sims/"
setwd(directory)
source("/Users/adwillis/research/alpha-aitchison/utility.R")

devtools::install_github("jacobbien/simulator")
library(simulator) 
library(dplyr)
library(devtools)
library(magrittr)
library(tidyverse)
library(dplyr)
library(foreach)
#  detach(package:plyr)
library(tidyr)
library(iNEXT)
library(grid)
library(gridExtra)

build("/Users/adwillis/software/breakaway/", vignettes=F)
install(pkg = "/Users/adwillis/software/breakaway/")
library(breakaway)

build("/Users/adwillis/research/alpha-aitchison/MicrobiomePack/")
install(pkg = "/Users/adwillis/research/alpha-aitchison/MicrobiomePack/")
library(MicrobiomePack)

build("/Users/adwillis/research/alpha-aitchison/Arbel/")
roxygen2::roxygenise("/Users/adwillis/research/alpha-aitchison/Arbel/")
install(pkg = "/Users/adwillis/research/alpha-aitchison/Arbel/")
library(Arbel)


install_github("zdk123/SpiecEasi")
library(SpiecEasi)

# Check 1: Sample covariance is a fine estimate as n gets large
# Test this by repeatedly simulating from model with fixed covariance
# and looking at distance between estimated covariance and truth as n gets large

## @knitr main
source("model_functions.R")
covariance_simulation <- new_simulation(name = "check-covariance",
                                        label = "Check that sample covariance converges to truth as n gets large") %>%
  generate_model(simulate_data_sigma, seed = 1234,
                 n = as.list(c(10, 28, 46, 62, 80, 98, 116)),
                 q = 5,
                 beta_sd = 1,
                 sigma_min = 0.01,
                 sigma_max = 10,
                 vary_along = c("n")) %>%
  simulate_from_model(nsim = 3)

source("method_functions.R")
covariance_simulation_run <- covariance_simulation %>% 
  run_method(list(get_covariance_rich, 
                  get_covariance_poor))

source("eval_functions.R")
covariance_simulation_evaluated <- covariance_simulation_run %>% 
  evaluate(list(mse_covariance, max_norm_covariance, mse_inverse_covariance, max_norm_inverse_covariance)) 

plot_eval_by(covariance_simulation_evaluated, "mse_covariance", varying = "n")
plot_eval_by(covariance_simulation_evaluated, "mse_inverse_covariance", varying = "n")
plot_eval_by(covariance_simulation_evaluated, "max_norm_covariance", varying = "n")
plot_eval_by(covariance_simulation_evaluated, "max_norm_inverse_covariance", varying = "n")
plot_eval_by(covariance_simulation_evaluated, "time", varying = "n")

q=5
A <- matrix(runif(q^2, min=-1, max=1), ncol=q) 
D <- diag(seq(from = 10, to=0.01, length.out=q))
sigma <- t(A) %*% D %*% A
mw20 <- make_w(mu=matrix(1, nrow = 20, ncol = 5), Sigma=sigma, mm = 1e5)
poor20 <- fit_aitchison(mw20, X=matrix(1, nrow = 20, ncol = 1), poorman=T)$sigma
rich20 <- fit_aitchison(mw20, X=matrix(1, nrow = 20, ncol = 1), poorman=F)$sigma
mw100 <- make_w(mu=matrix(1, nrow = 100, ncol = 5), Sigma=sigma, mm = 1e5)
poor100 <- fit_aitchison(mw100, X=matrix(1, nrow = 100, ncol = 1), poorman=T, in_parallel=TRUE,
                         EMiter = 6, 
                         EMburn = 4, 
                         MCiter = 200, 
                         MCburn = 100)$sigma
rich100 <- fit_aitchison(mw100, X=matrix(1, nrow = 100, ncol = 1), poorman=F, in_parallel=TRUE,
                         EMiter = 6, 
                         EMburn = 4, 
                         MCiter = 200, 
                         MCburn = 100)$sigma
poor20
rich20
poor100
rich100
sigma
#### They still seem to be identical
#### Maybe with larger sigma_max there will appear a difference?
### Ah, they are different on the inverse scale

###########################################################################
########## Now to use this to estimate variance ###########################
###########################################################################
source("make_lee_data.R")

source("model_functions.R")
covariance_simulation <- new_simulation(name = "check-alpha-variance",
                                        label = "Decide how to best get variance") %>%
  generate_model(simulate_data_sigma, seed = 1234,
                 q = as.list(seq(from = 5, to = 20, length.out = 4)),
                 beta_sd = 1,
                 vary_along = c("q")) %>%
  simulate_from_model(nsim = 100)
covariance_simulation


source("method_functions.R")
covariance_simulation_run <- covariance_simulation %>%
  run_method(c(sapply(c(0,1,2), make_parametric_variance_estimate),
               sapply(c(0,1,2), make_nonparametric_variance_estimate),
               sapply(c(0,1,2), make_delta_variance_estimate)))

source("eval_functions.R")
covariance_simulation_evaluated <- covariance_simulation_run %>% 
  evaluate(list(variance_difference_bc,
                variance_difference_euclidean,
                variance_difference_shannon_glassy, 
                variance_difference_shannon_altered,
                #variance_shannon_glassy,
                #variance_shannon_altered,
                variance_difference_simpson_glassy,
                variance_difference_simpson_altered)) 

ev_list <- simulator:::get_evals_list(covariance_simulation_evaluated)
evals_df <- rbind(as.data.frame(ev_list[1]),
                  as.data.frame(ev_list[2]))

qs <- evals_df$Model  %>% as.character %>% strsplit("_") %>% lapply(function(x) x[4])  %>% unlist %>% as.numeric
my_variance <- evals_df$Method %>% as.character %>% strsplit("_")%>% lapply(function(x) x[1])  %>% unlist 
my_covariance <- evals_df$Method  %>% as.character %>% strsplit("_")%>% lapply(function(x) x[2])  %>% unlist 

all_data_frame <- data.frame(evals_df, "q" = qs, "covariance" = my_covariance, "variance" = my_variance)


ggplot(all_data_frame,
       aes(x = q, y = variance_shannon_glassy, col = variance)) +
  geom_point() +
  facet_grid(.~covariance)

ggplot(all_data_frame,
       aes(x = q, y = variance_shannon_glassy, col = covariance, 
           pch = covariance)) +
  geom_point() +
  facet_grid(.~variance)


ggplot(data.frame(all_data_frame, "h" = paste(all_data_frame$q, all_data_frame$variance)),
       aes(x = h, y = variance_shannon_glassy, col = variance, 
           pch = covariance)) +
  geom_point() +
  geom_boxplot() +
  facet_grid(.~covariance)  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(data.frame(all_data_frame, "h" = paste(all_data_frame$q, all_data_frame$variance)) %>%
         filter(substr(variance, 1,2) %in% c("no", "pa")),
       aes(x = h, y = variance_difference_shannon_glassy, col = variance, 
           pch = covariance)) +
  geom_point() +
  geom_boxplot() +
  facet_grid(.~covariance) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



# pdf("shannon_variance_difference.pdf")
ggplot(data.frame(all_data_frame, "h" = paste(all_data_frame$q, all_data_frame$variance)),
       aes(x = h, y = variance_difference_shannon_glassy, col = variance, 
           pch = covariance)) +
  geom_point() +
  geom_boxplot() +
  facet_grid(.~covariance) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# dev.off()

ggplot(data.frame(all_data_frame, "h" = paste(all_data_frame$q, all_data_frame$variance)),
       aes(x = h, y = variance_difference_simpson_glassy, col = variance, 
           pch = covariance)) +
  geom_point() +
  geom_boxplot() +
  facet_grid(.~covariance) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


all_data_frame[, c("variance_shannon_glassy", "q", "variance", "covariance")] %>% 
  filter(q <= 10 & covariance == "stars") %>% 
  ggplot(aes(x = ))

all_data_frame %>% names

ggplot(data.frame(evals_df, qs), 
       aes(x=qs, y = variance_difference_shannon_glassy, col = Method)) +
  geom_point()

ggplot(data.frame(evals_df, qs), 
       aes(x=qs, y = variance_difference_shannon_altered, col = Method)) +
  geom_point()

ggplot(data.frame(evals_df, qs), 
       aes(x=qs, y = variance_shannon_altered, col = Method)) +
  geom_point()

ggplot(data.frame(evals_df, qs) %>% 
         filter(substr(Method, 1, 2) %in% c("no","pa") ), 
       aes(x=qs, y = variance_shannon_glassy, col = Method)) +
  geom_point() 

ggplot(data.frame(evals_df, qs), 
       aes(x=qs, y = variance_shannon_glassy, col = Method)) +
  geom_point() 
ggplot(data.frame(evals_df, qs), 
       aes(x=qs, y = time, col = Method)) +
  geom_point() 

plot_eval_by(covariance_simulation_evaluated, "variance_difference_shannon_altered", varying = "q")
plot_eval_by(covariance_simulation_evaluated, "variance_difference_shannon_glassy", varying = "q")
plot_eval_by(covariance_simulation_evaluated, "variance_shannon_altered", varying = "q")
plot_eval_by(covariance_simulation_evaluated, "variance_shannon_glassy", varying = "q")

#save.image("end180307.RData")




########## Start small
source("make_lee_data.R")
source("model_functions.R")
covariance_simulation_small <- new_simulation(name = "check-alpha-variance-smaller",
                                              label = "Decide how to best get variance (baby simulation)") %>%
  generate_model(simulate_data_sigma, seed = 1234,
                 q = as.list(seq(from = 20, to = 5, length.out = 2)),
                 beta_sd = 1,
                 vary_along = c("q")) %>%
  simulate_from_model(nsim = 2)#, index = 1
covariance_simulation_small

source("method_functions.R")
covariance_simulation_run_small <- covariance_simulation_small %>%
  run_method(c(sapply(c(0,1,2), make_parametric_variance_estimate),
               sapply(c(0,1,2), make_nonparametric_variance_estimate)))

source("eval_functions.R")
covariance_simulation_evaluated_small <- covariance_simulation_run_small %>% 
  evaluate(list(variance_difference_bc,
                variance_difference_euclidean,
                variance_difference_shannon_glassy, 
                variance_difference_shannon_altered,
                variance_difference_simpson_glassy,
                variance_difference_simpson_altered)) 

save_simulation(covariance_simulation_evaluated_small)

## add to it
source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")

covariance_simulation_evaluated_small %<>%
  simulate_from_model(nsim = 10, index = 2) %>%
  run_method(c(sapply(c(2,1,0), make_parametric_variance_estimate),
               sapply(c(2,1,0), make_nonparametric_variance_estimate)), 
             parallel = list(socket_names = 6))


source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")

covariance_simulation_evaluated_small %<>%
  simulate_from_model(nsim = 10, index = 2) %>%
  run_method(c(sapply(c(2,1,0), make_parametric_variance_estimate),
               sapply(c(2,1,0), make_nonparametric_variance_estimate)), 
             parallel = list(socket_names = 6))
covariance_simulation_evaluated_small


covariance_simulation_evaluated_small %<>%
  simulate_from_model(nsim = 10, index = 2) %>%
  run_method(c(sapply(c(2,1,0), make_parametric_variance_estimate),
               sapply(c(2,1,0), make_nonparametric_variance_estimate)), 
             parallel = list(socket_names = 6))
covariance_simulation_evaluated_small

covariance_simulation_evaluated_small %<>%
  simulate_from_model(nsim = 10, index = 3) %>%
  run_method(c(sapply(c(2,1,0), make_parametric_variance_estimate),
               sapply(c(2,1,0), make_nonparametric_variance_estimate)), 
             parallel = list(socket_names = 8))
covariance_simulation_evaluated_small

# 
# # Hint: The following code can be used to recreate the error, where 'met' is the method object
# # (i.e. met@name == "parametric_stars") and 'sim' is your simulation object:
# 
# sim = covariance_simulation_evaluated_small  
# m <- model(sim, subset = "aitchison-model/beta_sd_1/q_20")
# d <- draws(sim, subset = "aitchison-model/beta_sd_1/q_20", index = 1)
# d@draws
# .Random.seed <<- as.integer(c(407, -430235322, -2061546260, -149061901, -1614633419, -731355787, -1895136781))
# make_parametric_variance_estimate(0)@method(model = m, draw = d@draws$r1.1)
# make_parametric_variance_estimate(1)@method(model = m, draw = d@draws$r1.1)
# make_parametric_variance_estimate(2)@method(model = m, draw = d@draws$r1.1)
# make_parametric_variance_estimate(2)@method(model = m, draw = d@draws$r1.2)
# make_nonparametric_variance_estimate(0)@method(model = m, draw = d@draws$r1.1)
# 
# 

source("eval_functions.R")
covariance_simulation_evaluated_small <- covariance_simulation_evaluated_small %>% 
  evaluate(list(variance_difference_bc,
                variance_difference_euclidean,
                variance_difference_shannon_glassy, 
                variance_difference_shannon_altered,
                variance_difference_simpson_glassy,
                variance_difference_simpson_altered)) 


ev_list <- simulator:::get_evals_list(covariance_simulation_evaluated_small)
evals_df <- as.data.frame(ev_list)

evals_df$Run <- evals_df$Draw %>% as.character %>% strsplit("r") %>% lapply(function(x) x[2]) %>% lapply(substring, 1, 1) %>% unlist 
evals_df <- evals_df[evals_df$Run == "2", ]   
qs <- evals_df$Model  %>% as.character %>% strsplit("_") %>% lapply(function(x) x[4])  %>% unlist %>% as.numeric
my_variance <- evals_df$Method %>% as.character %>% strsplit("_")%>% lapply(function(x) x[1])  %>% unlist 
my_covariance <- evals_df$Method  %>% as.character %>% strsplit("_")%>% lapply(function(x) x[2])  %>% unlist 

all_data_frame <- data.frame(evals_df, "q" = qs, "covariance" = my_covariance, "variance" = my_variance)
ordered <- data.frame(all_data_frame, "h" = paste(all_data_frame$q, all_data_frame$variance, sep = "_") %>% as.character)
ordered <- ordered[order(ordered$h, decreasing=T), ]

ggplot(ordered,
       aes(x = h, y = variance_difference_simpson_glassy, col = variance, 
           pch = covariance)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(.~covariance) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))





bc <- ggplot(ordered,
             aes(x = h, y = variance_difference_bc, col = variance, 
                 pch = covariance)) +
  # geom_boxplot() +
  geom_violin()+
  facet_grid(.~covariance) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

euc <- ggplot(ordered,
              aes(x = h, y = variance_difference_euclidean, col = variance, 
                  pch = covariance)) +
  # geom_boxplot() +
  geom_violin() +
  facet_grid(.~covariance) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

grid.arrange(bc, euc)




########## ########## ########## 
########## Get bigger ########## 
########## ########## ########## 

source("make_lee_data.R")
source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")

sim <- new_simulation(name = "check-alpha-variance-bigger",
                      label = "Decide how to best get variance (big simulation)") %>%
  generate_model(simulate_data_sigma, seed = 123,
                 q = as.list(seq(from = 20, to = 5, length.out = 4)),
                 beta_sd = 1,
                 vary_along = c("q")) %>%
  simulate_from_model(nsim = 4, index = 1:5) # 20 draws from each model

sim %<>% run_method(c(sapply(c(0,1,2), make_parametric_variance_estimate),
                      sapply(c(0,1,2), make_nonparametric_variance_estimate)), 
                    parallel = list(socket_names = 5)) 
# just ran above line

sim %<>% evaluate(list(variance_difference_bc,
                       variance_difference_euclidean,
                       variance_difference_shannon_glassy, 
                       variance_difference_shannon_altered,
                       variance_difference_simpson_glassy,
                       variance_difference_simpson_altered)) 

save_simulation(sim)

## add another 60 draws to it
sim %<>%
  simulate_from_model(nsim = 7, index = 6:10) %>%
  run_method(c(sapply(c(0,1,2), make_parametric_variance_estimate),
               sapply(c(0,1,2), make_nonparametric_variance_estimate)), 
             parallel = list(socket_names = 5))
sim %>% 
  evaluate(list(variance_difference_bc,
                variance_difference_euclidean,
                variance_difference_shannon_glassy, 
                variance_difference_shannon_altered,
                variance_difference_simpson_glassy,
                variance_difference_simpson_altered))

## add another 20 draws to it
sim %<>%
  simulate_from_model(nsim = 4, index = 11:15) %>%
  run_method(c(sapply(c(0,1,2), make_parametric_variance_estimate),
               sapply(c(0,1,2), make_nonparametric_variance_estimate)), 
             parallel = list(socket_names = 5)) 

sim %<>%  evaluate(list(variance_difference_bc,
                        variance_difference_euclidean,
                        variance_difference_shannon_glassy, 
                        variance_difference_shannon_altered,
                        variance_difference_simpson_glassy,
                        variance_difference_simpson_altered, 
                        variance_shannon_glassy,
                        variance_bc))

save_simulation(sim)
ev_list <- simulator:::get_evals_list(sim)
evals_df <- as.data.frame(ev_list)

# evals_df$Run <- evals_df$Draw %>% as.character %>% strsplit("r") %>% lapply(function(x) x[2]) %>% lapply(substring, 1, 1) %>% unlist 
# evals_df <- evals_df[evals_df$Run == "2", ]   
qs <- evals_df$Model  %>% as.character %>% strsplit("_") %>% lapply(function(x) x[4])  %>% unlist %>% as.numeric
my_variance <- evals_df$Method %>% as.character %>% strsplit("_")%>% lapply(function(x) x[1])  %>% unlist 
my_covariance <- evals_df$Method  %>% as.character %>% strsplit("_")%>% lapply(function(x) x[2])  %>% unlist 

all_data_frame <- data.frame(evals_df, "q" = qs, "covariance" = my_covariance, "variance" = my_variance)
my_h <- paste(all_data_frame$q, all_data_frame$variance, sep = "_") %>% as.factor
print(levels(my_h))
my_h = factor(my_h, levels(my_h)[c(8, 7, 2, 1, 4, 3,6,5)])
print(levels(my_h))

ordered <- data.frame(all_data_frame, "h" = my_h)
# ordered <- ordered[order(ordered$h, decreasing=F), ]

### 


simpson_glassy <- ggplot(ordered,
                         aes(x = h, y = variance_difference_simpson_glassy, col = variance, 
                             pch = covariance)) +
  geom_violin() + facet_grid(.~covariance) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
simpson_glassy

base_plot <- ggplot(ordered,
                    aes(x = h, col = variance, 
                        pch = covariance)) +
  facet_grid(.~covariance, labeller = as_labeller(c("diagonal" = "Diagonal\ncovariance",
                                                    "naive" = "Unstructured\ncovariance",
                                                    "stars" = "Graphical\nlasso"))) + 
  scale_color_manual(name = "Variance estimation method", labels = c("Nonparametric", "Parametric"),
                     values = c("orange", "darkblue")) +
  xlab("Number of taxa") + 
  scale_x_discrete(labels = rep(x=c("   5", "", "   10", "", "   15", "", "   20", ""), 1)) 

shannon_glassy_difference <- base_plot + 
  geom_violin(aes(x = h, variance_difference_shannon_glassy)) + 
  ylab("Difference between\nestimated variance and\ntrue variance (Shannon)") 
shannon_glassy_difference

shannon_glassy <- base_plot +
  geom_violin(aes(x = h, y = variance_shannon_glassy)) + 
  ylab("Estimated variance\n(Shannon)") 
shannon_glassy

bc <-  base_plot +
  geom_violin( aes(x = h, y = variance_bc)) + 
  ylab("Estimated variance\n(Bray-Curtis)")
bc

bc_difference <- base_plot +
  geom_violin( aes(x = h, y = variance_difference_bc)) + 
  ylab("Difference between\nestimated variance and\ntrue variance (Bray-Curtis)") 
bc_difference

ggpubr::ggarrange(shannon_glassy_difference, 
                  shannon_glassy, 
                  bc_difference, bc, 
                  ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

ggsave("variance_estimation.pdf", width=8, height = 5)

##### 

ggplot(ordered,
       aes(x = Method, col = variance, 
           pch = covariance)) +
  geom_violin(aes(x = Method, y = variance_shannon_glassy)) + 
  facet_grid(~q)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))


# load("end180313.RData")

#save.image("end180313.RData")
ordered$h %>% table


# sim0 = covariance_simulation_evaluated_small
# m <- model(sim0, subset = "aitchison-model/beta_sd_1/q_20")
# d <- draws(sim0, subset = "aitchison-model/beta_sd_1/q_20", index = 1)
# d@draws
# .Random.seed <<- as.integer(c(407, -430235322, -2061546260, -149061901, -1614633419, -731355787, -1895136781))
# arbel@method(model = m, draw = d@draws$r1.1)


###### 


source("make_lee_data.R")
source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")

sim_60 <- new_simulation(name = "check-alpha-variance-bigger-again",
                         label = "Decide how to best get variance (bigger simulation)") %>%
  generate_model(simulate_data_sigma_rich_structure, seed = 123,
                 q = as.list(c(60,40,20)),
                 beta_sd = 5,
                 vary_along = c("q")) %>%
  simulate_from_model(nsim = 2, index = 1:5) # 10 draws from each model
sim_60 %<>% run_method(c(sapply(c(0,1,2), make_parametric_variance_estimate),
                         sapply(c(0,1,2), make_nonparametric_variance_estimate)), 
                       parallel = list(socket_names = 5)) 
system('say "done"')
save_simulation(sim_60)
sim_60 %<>%  evaluate(list(variance_difference_bc,
                           variance_difference_euclidean,
                           variance_difference_shannon_glassy, 
                           variance_difference_shannon_altered,
                           variance_difference_simpson_glassy,
                           variance_difference_simpson_altered, 
                           variance_shannon_glassy,
                           variance_bc))


sim_60 <- sim_60 %>%
  simulate_from_model(nsim = 18, index = 6:10) %>%
  run_method(c(sapply(c(0,1,2), make_parametric_variance_estimate),
               sapply(c(0,1,2), make_nonparametric_variance_estimate)), 
             parallel = list(socket_names = 5))

sim_60 %<>%  evaluate(list(variance_difference_bc,
                           variance_difference_euclidean,
                           variance_difference_shannon_glassy, 
                           variance_difference_shannon_altered,
                           variance_difference_simpson_glassy,
                           variance_difference_simpson_altered, 
                           variance_shannon_glassy,
                           variance_bc))
save_simulation(sim_60)

ev_list <- simulator:::get_evals_list(sim_60)
evals_df <- as.data.frame(ev_list)
qs <- evals_df$Model  %>% as.character %>% strsplit("_") %>% lapply(function(x) x[4])  %>% unlist %>% as.numeric
my_variance <- evals_df$Method %>% as.character %>% strsplit("_")%>% lapply(function(x) x[1])  %>% unlist 
my_variance[my_variance == "nonparametric"] <- "np"
my_variance[my_variance == "parametric"] <- "P"
my_covariance <- evals_df$Method  %>% as.character %>% strsplit("_")%>% lapply(function(x) x[2])  %>% unlist 
all_data_frame <- data.frame(evals_df, "q" = qs, "covariance" = my_covariance, "variance" = my_variance)
ordered <- data.frame(all_data_frame)

#### variance figure
base_variance <- ggplot(ordered,
                        aes(x = Method, col = variance, fill = variance, 
                            alpha = covariance)) +
  facet_grid(~q, labeller = as_labeller(c("20" = "20 taxa",
                                          "40" = "40 taxa",
                                          "60" = "60 taxa")))+ 
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # theme_bw() +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(name = "Variance: ", labels = c("Nonparametric", "Parametric"),
                     values = c("orange", "darkblue")) +
  scale_fill_manual(name = "Variance: ", labels = c("Nonparametric", "Parametric"),
                    values = c("orange", "darkblue")) +
  scale_alpha_manual(guide = F,
                     values = c(0.3, 0.5, 0.7)) + xlab("") +
  scale_x_discrete(labels = c("Diagonal", "Naive", "G. lasso", "Diagonal", "Naive", "G. lasso")) 
var_shannon <- base_variance + geom_boxplot(aes(x = Method, y = variance_shannon_glassy)) ; var_shannon

var_shannon <- base_variance + ylab("Estimated variance\n(Shannon)") +
  geom_boxplot(aes(x = Method, y = variance_shannon_glassy)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.3))

var_bc <- base_variance + ylab("Estimated variance\n(Bray-Curtis)") +
  geom_boxplot(aes(x = Method, y = variance_bc))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.15))

var_diff_shannon <- base_variance + ylab("Difference between\nestimated variance and\ntrue variance (Shannon)") + 
  geom_hline(yintercept=0, col = "grey") + 
  geom_boxplot(aes(x = Method, y = variance_difference_shannon_glassy))+
  scale_y_continuous(expand = c(0, 0), limits = c(-2.2,2.2))

var_diff_bc <- base_variance + ylab("Difference between\nestimated variance and\ntrue variance (Bray-Curtis)") +
  geom_hline(yintercept=0, col = "grey") + 
  geom_boxplot(aes(x = Method, y = variance_difference_bc)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-.15,0.15)) 
var_diff_bc


ggpubr::ggarrange(var_shannon, 
                  var_diff_shannon, 
                  var_bc, var_diff_bc, 
                  ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

ggsave("larger_variance_estimation_clean.pdf", width=8, height = 6)
save.image("end_variance_plots.RData")
