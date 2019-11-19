# Amy Willis, 19 November 2019
# A script to create Figure 4 of the DivNet paper
# sigma_max vs MSE

library(simulator) # this file was created under simulator version 0.2.0
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
source("../simulator_files/method_functions.R")
source("../simulator_files/method_functions_secondary.R")
source("../simulator_files/eval_functions.R")


###########################################################################
########################## alpha diversity -- sigma_max ###############
###########################################################################
sim_alpha_sigma_max <- new_simulation(name = "sigma-max",
                                      label = "Estimating diversity varying sigma_max") %>%
  generate_model(make_aitchison, 
                 seed = 191119,
                 n = 20,
                 p = 2, 
                 q = 20,
                 sigma_min = 0.01,
                 sigma_max = as.list(c(0.1, 5, 10, 15)),
                 beta_sd = 1,
                 vary_along = c("sigma_max")) %>%
  simulate_from_model(nsim = 100)  

sim_alpha_sigma_max_methods <- sim_alpha_sigma_max %>%
  run_method(list(chao_shen,
                  divnet_method,
                  arbel,
                  inext,
                  plug,
                  plug_add_half), 
             parallel = list(socket_names = 6)) 
saveRDS(sim_alpha_sigma_max_methods, file="sim_alpha_sigma_max_methods.RDS")

sim_alpha_sigma_max_evals <- sim_alpha_sigma_max_methods %>% 
  evaluate(list(mse_loss_shannon, 
                mse_loss_simpson, 
                mse_loss_euclidean, 
                mse_loss_bray_curtis, 
                mse_loss_bray_curtis_subset, 
                mse_loss_euclidean_subset)) 

my_data_frame <- sim_alpha_sigma_max_evals %>% evals %>% as.data.frame

sms <- my_data_frame$Model %>% 
  as.character %>% 
  strsplit("/") %>% 
  lapply(function(x) x[6]) %>% 
  unlist %>% 
  strsplit(split="_") %>% 
  lapply(function(x) x[3]) %>% 
  unlist %>% 
  as.numeric

full_data_frame <- data.frame(my_data_frame, 
                              "sigma_max" = sms, 
                              "sigma_max_char" = as.factor(sms),
                              "sigma_max_method" = paste(sms, "_", my_data_frame$Method, sep = ""))

full_data_frame$Method %>% unique
full_data_frame$Method <- plyr::revalue(full_data_frame$Method, 
                                        c("divnet"="Proposed",
                                          "arbel" = "Arbel et. al",
                                          "chao-shen" = "Chao & Shen",
                                          "inext" = "iNEXT",
                                          "plug-in" = "Multinomial MLE",
                                          "plug-in-add-half" = "Zero-replace"))

full_data_frame$Method %<>% relevel("Multinomial MLE")
full_data_frame$Method %<>% relevel("Zero-replace")
full_data_frame$Method %<>% relevel("Proposed")

full_data_frame %<>%
  arrange(sigma_max)

full_data_frame %>%
  select(Method, mse_loss_shannon, sigma_max_char, sigma_max) %>% 
  as_tibble %>%
  filter(Method == "Arbel et. al", sigma_max > 1)
shannon_plot <- ggplot(full_data_frame, 
                       aes(x = sigma_max_char, 
                           y = mse_loss_shannon, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab(label=expression(sigma[max]~": Strength of cooccurrences")) +
  ylab("MSE: Shannon index")  +
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00", "#56B4E9",  "black"))  
shannon_plot

simpson_plot <- ggplot(full_data_frame, 
                       aes(x = sigma_max_char, y = mse_loss_simpson, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab(label=expression(sigma[max]~": Strength of cooccurrences")) +
  ylab("MSE: Simpson index") +
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00", "black")) 
simpson_plot

bc_plot <- ggplot(full_data_frame, 
                  aes(x = sigma_max_char, y = mse_loss_bray_curtis, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab(label=expression(sigma[max]~": Strength of cooccurrences")) +
  ylab("Log MSE: Bray-Curtis index") +
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00", "#56B4E9"))  
bc_plot

euc_plot <- ggplot(full_data_frame, 
                   aes(x = sigma_max_char, 
                       y = mse_loss_euclidean_subset, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab(label=expression(sigma[max]~": Strength of cooccurrences")) +
  ylab("Log MSE: Euclidean index") +
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00", "#56B4E9")) 
euc_plot

time_plot <- ggplot(full_data_frame, 
                    aes(x = sigma_max_char, y = time, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab(label=expression(sigma[max]~": Strength of cooccurrences")) +
  ylab("Mean computational time (seconds)")+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00", "#56B4E9", "black")) 
time_plot

ggpubr::ggarrange(shannon_plot, 
                  simpson_plot, 
                  bc_plot, 
                  euc_plot, ncol=2, nrow=2, 
                  common.legend = T, legend="right")

