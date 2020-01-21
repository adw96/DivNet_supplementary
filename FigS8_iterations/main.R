# Amy Willis, 19 November 2019
# A script to create Figure 4 of the DivNet paper
# n vs MSE

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


#### vary sigma_max

source("../simulator_files/model_functions.R")
source("../simulator_files/eval_functions.R")
source("../simulator_files/method_functions_iterations.R")

sim_runtime <- new_simulation(name = "divnet-n",
                              label = "Estimating diversity varying n") %>%
  generate_model(make_aitchison, seed = 200113,
                 n = 10,
                 p = 2, 
                 q = 20,
                 sigma_min = 0.01,
                 sigma_max = 5,
                 beta_sd = 1) %>%
  simulate_from_model(nsim = 10) 

sim_runtime <- sim_runtime %>%
  run_method(list(arbel_500,
                  arbel_100,
                  arbel_10,
                  divnet_50,
                  divnet_20,
                  divnet_6)) %>% 
  evaluate(list(mse_loss_shannon, 
                estimated_shannon, 
                true_shannon)) 

my_data_frame <- sim_runtime %>% evals %>% as.data.frame
my_data_frame

my_data_frame %>%
  ggplot(aes(x = Method, 
             y = estimated_shannon, col = Method)) +
  geom_boxplot(outlier.size=0.1)

my_data_frame %>%
  ggplot(aes(x = Method, 
             y = mse_loss_shannon, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab(label="Method") +
  ylab("MSE: Shannon index")  +
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) 


####

source("../simulator_files/model_functions.R")
source("../simulator_files/eval_functions.R")
source("../simulator_files/method_functions_iterations.R")
sim_runtime_all_arbel <- new_simulation(name = "divnet-n",
                                        label = "Estimating diversity varying n") %>%
  generate_model(make_aitchison, seed = 200113,
                 n = 40,
                 p = 2, 
                 q = 100,
                 sigma_min = 0.01,
                 sigma_max = 5,
                 beta_sd = 1) %>%
  simulate_from_model(nsim = 10) %>%
  run_method(list(divnet_10,
                  arbel_2000,
                  arbel_1000,
                  arbel_500,
                  arbel_250,
                  divnet_6,
                  divnet_20)) %>% 
  evaluate(list(mse_loss_shannon, 
                estimated_shannon, 
                true_shannon)) 
# saveRDS(sim_runtime_all_arbel, "sim_runtime_all_arbel.RDS")

# write_csv(sim_runtime_all_arbel %>% 
#   evals %>% 
#   as.data.frame, "sim_runtime_all_arbel_evals_df.csv")
sim_runtime_all_arbel_evals_df <- read_csv("sim_runtime_all_arbel_evals_df.csv")

cleaned_df <- sim_runtime_all_arbel_evals_df %>% 
  separate(Method, into=c("Method", "Iterations"), sep="(?<=[A-Za-z])(?=[0-9])") %>%
  mutate(Method = str_replace(Method, "divnet", "DivNet with")) %>%
  mutate(Method = str_replace(Method, "arbel", "Arbel et al with")) %>%
  mutate(Iterations = paste(Iterations, 
                            ifelse(Method == "DivNet with", 
                                   "EM iterations", "MC iterations"))) %>%
  unite("Method", Method:Iterations, sep = " ") %>%
  mutate(Method = factor(Method, levels=c("DivNet with 6 EM iterations",
                                          "DivNet with 10 EM iterations",
                                          "DivNet with 20 EM iterations",
                                          "Arbel et al with 250 MC iterations",
                                          "Arbel et al with 500 MC iterations",
                                          "Arbel et al with 1000 MC iterations",
                                          "Arbel et al with 2000 MC iterations")))
mse_plot <- cleaned_df %>%
  ggplot(aes(x = Method, 
             y = mse_loss_shannon, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab(label="Method") +
  ylab("MSE: Shannon index")  +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) 

time_plot <- cleaned_df %>% 
  ggplot(aes(x = Method, 
             y = time, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab(label="Method") +
  ylab("Computation time (seconds)")  +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) 
ggpubr::ggarrange(mse_plot, time_plot, ncol = 2, common.legend=TRUE, legend="right")
# ggsave("time-mse-comparison.pdf", width=8, height = 3)

cleaned_df %>%
  group_by(Method) %>%
  filter(Method %in% c("DivNet with 6 EM iterations",
                       "Arbel et al with 500 MC iterations")) %>%
  summarise(median(time))

### just to check

source("../simulator_files/model_functions.R")
source("../simulator_files/eval_functions.R")
source("../simulator_files/method_functions_iterations.R")
mini_div <- new_simulation(name = "divnet-n",
                                        label = "Estimating diversity varying n") %>%
  generate_model(make_aitchison, seed = 200113,
                 n = 10,
                 p = 2, 
                 q = 10,
                 sigma_min = 0.01,
                 sigma_max = 5,
                 beta_sd = 1) %>%
  simulate_from_model(nsim = 10) %>%
  run_method(list(#divnet_20
                  #divnet_10,
                  divnet_2,
                  divnet_6
                  )) %>% 
  evaluate(list(mse_loss_shannon, 
                estimated_shannon, 
                true_shannon)) 

mini_div %>% 
    evals %>%
    as.data.frame %>%
  ggplot(aes(x = Method, 
             y = mse_loss_shannon, col = Method)) +
  geom_violin(outlier.size=0.1)

mini_div %>% 
  evals %>%
  as_tibble %>%
  pivot_wider(id_cols=Draw, names_from=Method, values_from=mse_loss_shannon) %>%
  ggplot(aes(x = divnet2, 
             y = divnet6)) + geom_point() + geom_abline(aes(slope=1, intercept=0))+
  scale_x_log10() + scale_y_log10()
