library(devtools)
library(simulator) 
library(dplyr)
library(devtools)
library(magrittr)
library(tidyverse)
library(dplyr)
library(foreach)
library(tidyr)
library(iNEXT)
library(grid)
library(gridExtra)
library(breakaway)
library(ggplot2)

install_github("zdk123/SpiecEasi")
library(SpiecEasi)

devtools::build("/Users/adwillis/software/DivNet/", vignettes=F)
devtools::install("/Users/adwillis/software/DivNet/")
library(DivNet)

###### Resubmission


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
  # simulate_from_model(nsim = 2, index = 1:5) # 10 draws from each model
  simulate_from_model(nsim = 2, index = 1) # 10 draws from each model

sim_60 %<>% run_method(sapply(c(0,1), make_nonparametric_variance_estimate))


source("method_functions.R")
make_nonparametric_variance_estimate(0)@method(model = m, draw = d@draws$r1.1, cov_method = "default")
make_nonparametric_variance_estimate(1)@method(model = m, draw = d@draws$r1.1, cov_method = "diagonal")
make_nonparametric_variance_estimate(2)@method(model = m, draw = d@draws$r1.1, cov_method = "stars")


sim_60 %>% run_method(sapply(c(0,1,2), make_parametric_variance_estimate))

sim_60 %<>% run_method(c(sapply(c(0,1,2), make_parametric_variance_estimate),
                         sapply(c(0,1,2), make_nonparametric_variance_estimate)), 
                       parallel = list(socket_names = 5)) 
system('say "done"')
save_simulation(sim_60)

sim_60 %<>%  evaluate(list(variance_difference_bc,
                           variance_difference_shannon_glassy, 
                           variance_difference_shannon_altered,
                           variance_shannon_glassy,
                           variance_bc))


sim_60 <- sim_60 %>%
  simulate_from_model(nsim = 18, index = 6:10) %>%
  run_method(c(sapply(c(0,1,2), make_parametric_variance_estimate),
               sapply(c(0,1,2), make_nonparametric_variance_estimate)), 
             parallel = list(socket_names = 5))

sim_60 %<>%  evaluate(list(variance_difference_bc,
                           variance_difference_shannon_glassy, 
                           variance_difference_shannon_altered,
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
