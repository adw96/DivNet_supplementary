# Amy Willis, 2 Dec 2019
# A script to create Figure 9 of the DivNet paper
# 

library(simulator) 
library(devtools)
library(magrittr)
library(tidyverse)
library(foreach)
library(iNEXT)
library(grid)
library(gridExtra)

install_github("adw96/breakaway")
library(breakaway)

install_github("adw96/DivNet")
library(DivNet)

install_github("jarbel/Dep-GEM")
library(DepGEM)

source("../simulator_files/model_functions_nonlinear.R")
source("../simulator_files/method_functions_nonlinear.R")
source("../simulator_files/eval_functions_nonlinear.R")

nonlinear_simulation_1 <- new_simulation(name = "nonlinear-simulation",
                                       label = "Evaluating robustness to misspecification") %>%
  generate_model(make_aitchison_nonlinear_function1, 
                 seed = 191202,
                 n = 20,
                 q = 40,
                 gamma1 = as.list(c(0.125, 0.25, 0.5)),
                 sigma_min = 0.01,
                 sigma_max = 5,
                 vary_along = c("gamma1")) %>%
  simulate_from_model(nsim = 100) 
nonlinear_simulation_1_methods <- nonlinear_simulation_1 %>%
  run_method(list(divnet_method_nonlinear, 
                  divnet_method_linear, 
                  inext,
                  chao_shen,
                  plug, 
                  arbel_linear,
                  arbel_quad
                  )) 
nonlinear_simulation_1_evals <- nonlinear_simulation_1_methods %>% 
  evaluate(list(mse_loss_shannon)) 

nonlinear_results <- right_join(nonlinear_simulation_1_methods %>% model %>% as.data.frame, 
                                nonlinear_simulation_1_evals %>% evals %>% as.data.frame, 
                                by = c("name" = "Model"))

sim1 <- nonlinear_results %>%
  mutate(gamma1 = as.character(gamma1)) %>%
  mutate(Method = plyr::revalue(Method, 
                                c("divnet-lin"="Proposed (Linear)",
                                  "divnet-nonlin"="Proposed (Quadratic)",
                                  "arbel-linear" = "Arbel et. al (Linear)",
                                  "arbel-quad" = "Arbel et. al (Quadratic)",
                                  "chao-shen" = "Chao & Shen",
                                  "inext" = "iNEXT",
                                  "plug-in" = "Multinomial MLE"))) %>%
  ggplot(aes(x = gamma1, y = mse_loss_shannon, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab(label=expression(gamma)) +
  ylab("MSE: Shannon") +
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  scale_color_manual(breaks = c("Proposed (Linear)", 
                                "Proposed (Quadratic)", 
                                "Multinomial MLE",
                                "Arbel et. al (Linear)", 
                                "Arbel et. al (Quadratic)",
                                "Chao & Shen", 
                                "iNEXT"), 
                     values = c("Proposed (Linear)" = "purple",
                                "Proposed (Quadratic)" = "blue",
                                "Arbel et. al (Linear)" = "#999000", 
                                "Arbel et. al (Quadratic)" = "#999999", 
                                "Chao & Shen" = "#E69F00",
                                "iNEXT" = "#56B4E9",
                                "Multinomial MLE" = "red")) 


sim1



##### Simulation 2
source("../simulator_files/model_functions_nonlinear.R")
source("../simulator_files/method_functions_nonlinear.R")
source("../simulator_files/eval_functions_nonlinear.R")

nonlinear_simulation_2 <- new_simulation(name = "nonlinear-simulation2",
                                         label = "Evaluating robustness to misspecification 2") %>%
  generate_model(make_aitchison_nonlinear_function2, 
                 seed = 191202,
                 n = 20,
                 q = 40,
                 gamma1 = as.list(c(0.125, 0.25, 0.32)),
                 sigma_min = 0.01,
                 sigma_max = 5,
                 vary_along = c("gamma1")) %>%
  simulate_from_model(nsim = 100) 
nonlinear_simulation_2_methods <- nonlinear_simulation_2 %>%
  run_method(list(divnet_method_nonlinear, 
                  divnet_method_linear, 
                  inext,
                  chao_shen,
                  plug, 
                  arbel_linear,
                  arbel_quad
  )) 
nonlinear_simulation_2_evals <- nonlinear_simulation_2_methods %>% 
  evaluate(list(mse_loss_shannon)) 

nonlinear_results_2 <- right_join(nonlinear_simulation_2_methods %>% model %>% as.data.frame, 
                                nonlinear_simulation_2_evals %>% evals %>% as.data.frame, 
                                by = c("name" = "Model"))

# write_csv(nonlinear_results, path="nonlinear_results_1.csv")
# write_csv(nonlinear_results_2, path="nonlinear_results_2.csv")
nonlinear_results_2$Method %>% unique

sim2 <- nonlinear_results_2 %>%
  mutate(gamma1 = as.character(gamma1)) %>%
  mutate(Method = plyr::revalue(Method, 
                                c("divnet-lin"="Proposed (Linear)",
                                  "divnet-nonlin"="Proposed (Quadratic)",
                                  "arbel-linear" = "Arbel et. al (Linear)",
                                  "arbel-quad" = "Arbel et. al (Quadratic)",
                                  "chao-shen" = "Chao & Shen",
                                  "inext" = "iNEXT",
                                  "plug-in" = "Multinomial MLE"))) %>%
  ggplot(aes(x = gamma1, y = mse_loss_shannon, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab(label=expression(gamma)) +
  ylab("MSE: Shannon") +
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  scale_color_manual(breaks = c("Proposed (Linear)", 
                                "Proposed (Quadratic)", 
                                "Multinomial MLE",
                                "Arbel et. al (Linear)", 
                                "Arbel et. al (Quadratic)",
                                "Chao & Shen", 
                                "iNEXT"), 
                     values = c("Proposed (Linear)" = "purple",
                                "Proposed (Quadratic)" = "blue",
                                "Arbel et. al (Linear)" = "#999000", 
                                "Arbel et. al (Quadratic)" = "#999999", 
                                "Chao & Shen" = "#E69F00",
                                "iNEXT" = "#56B4E9",
                                "Multinomial MLE" = "red")) 


sim2

# source("decide-true-model.R")
# ggpubr::ggarrange(quad_y, quad_relative, exp_quad, sim1, 
#                   exp_y, exp_relative, exp_shannon, sim2, nrow = 2, ncol = 4, 
#                   # common.legend=TRUE, legend="right",
#                   widths=c(1,1,1,2,1,1,1,2))
# ggsave("nonlinear-fig-n-100.pdf", width = 15, height = 5, units="in")

ggpubr::ggarrange(sim1 + ggtitle("Quadratic trend"), 
                  sim2 + ggtitle("Exponential trend"), nrow = 2,
                  common.legend=TRUE, legend="right")
ggsave("nonlinear-fig-n-100-results.pdf", width = 8, height = 5, units="in")

# ggpubr::ggarrange(quad_y, quad_relative, exp_quad, sim1, 
#                   exp_y, exp_relative, exp_shannon, sim2, nrow = 2, ncol = 4, 
#                   widths=c(1,1,1,2,1,1,1,2))
# ggsave("nonlinear-fig-n-10.pdf", width = 15, height = 5, units="in")
