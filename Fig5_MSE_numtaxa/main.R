# Amy Willis, 19 November 2019
# A script to create Figure 4 of the DivNet paper
# q vs MSE

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


#### vary q

source("../simulator_files/model_functions.R")
source("../simulator_files/method_functions.R")
source("../simulator_files/method_functions_secondary.R")
source("../simulator_files/eval_functions.R")


###########################################################################
########################## alpha diversity -- q ###############
###########################################################################
sim_alpha_q <- new_simulation(name = "sim-q",
                              label = "Estimating diversity varying q") %>%
  generate_model(make_aitchison, seed = 1234567,
                 n = 20,
                 p = 2, 
                 q = as.list(c(100, 60, 20)),
                 sigma_min = 0.01,
                 sigma_max = 5,
                 beta_sd = 1,
                 vary_along = c("q")) %>%
  simulate_from_model(nsim = 100) %>%
  run_method(list(chao_shen,
                  divnet_method,
                  arbel,
                  inext,
                  plug,
                  plug_add_half), 
             parallel = list(socket_names = 6)) %>% 
  evaluate(list(mse_loss_shannon, 
                mse_loss_simpson, 
                mse_loss_euclidean, 
                mse_loss_bray_curtis, 
                mse_loss_bray_curtis_subset, 
                mse_loss_euclidean_subset)) 


my_data_frame <- right_join(sim_alpha_q %>% model %>% as.data.frame, 
                            sim_alpha_q %>% evals %>% as.data.frame, 
                            by = c("name" = "Model"))
full_data_frame <- my_data_frame %>%
  as_tibble %>%
  separate(col=name, sep="/", into = c("model", "beta", "n", "p", 
                                       "q", "sigma_max", "sigma_min")) %>%
  separate(col = q, sep = "_", into = c("n1", "q_char")) %>%
  select(-n1) %>%
  mutate(Method = plyr::revalue(Method, 
                                c("divnet"="Proposed",
                                  "arbel" = "Arbel et. al",
                                  "chao-shen" = "Chao & Shen",
                                  "inext" = "iNEXT",
                                  "plug-in" = "Multinomial MLE",
                                  "plug-in-add-half" = "Zero-replace")))%>%
  mutate(Method = relevel(Method, "Multinomial MLE")) %>%
  mutate(Method = relevel(Method, "Zero-replace")) %>%
  mutate(Method = relevel(Method, "Proposed"))

# write_csv(full_data_frame, path="q_df.csv")
full_data_frame <- read_csv("q_df.csv")
full_data_frame %<>% 
  mutate(Method = as.factor(Method)) %>%
  # mutate(Method = plyr::revalue(Method, 
  #                               c("divnet"="Proposed",
  #                                 "arbel" = "Arbel et. al",
  #                                 "chao-shen" = "Chao & Shen",
  #                                 "inext" = "iNEXT",
  #                                 "plug-in" = "Multinomial MLE",
  #                                 "plug-in-add-half" = "Zero-replace")))%>%
  mutate(Method = relevel(Method, "Multinomial MLE")) %>%
  mutate(Method = relevel(Method, "Zero-replace")) %>%
  mutate(Method = relevel(Method, "Proposed")) %>%
  mutate(q_char = factor(q_char, levels=c("20", "60", "100"))) 

shannon_plot <- ggplot(full_data_frame, 
                       aes(x = q_char, 
                           y = mse_loss_shannon, 
                           col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab("Q: Number of taxa") +
  ylab("MSE: Shannon index")  +
  theme_bw() + 
  scale_color_manual(breaks = c("Proposed", 
                                "Multinomial MLE",
                                "Arbel et. al", 
                                "Chao & Shen", 
                                "iNEXT",
                                "Zero-replace"), 
                     values = c("Proposed" = "blue",
                                "Arbel et. al" = "#999999", 
                                "Chao & Shen" = "#E69F00",
                                "iNEXT" = "#56B4E9",
                                "Multinomial MLE" = "red",
                                "Zero-replace" = "black"))  +
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) 
shannon_plot

simpson_plot <- ggplot(full_data_frame, 
                       aes(x = q_char, 
                           y = mse_loss_simpson, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab("Q: Number of taxa") +
  ylab("MSE: Simpson index")  +
  theme_bw() + 
  scale_color_manual(breaks = c("Proposed", 
                                "Multinomial MLE",
                                "Arbel et. al", 
                                "Chao & Shen", 
                                "iNEXT",
                                "Zero-replace"), 
                     values = c("Proposed" = "blue",
                                "Arbel et. al" = "#999999", 
                                "Chao & Shen" = "#E69F00",
                                "iNEXT" = "#56B4E9",
                                "Multinomial MLE" = "red",
                                "Zero-replace" = "black"))  +
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) 
simpson_plot

bc_plot <- ggplot(full_data_frame, 
                  aes(x = q_char, 
                      y = mse_loss_bray_curtis, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab("Q: Number of taxa") + 
  ylab("MSE: Bray-Curtis index") +
  theme_bw() + 
  scale_color_manual(breaks = c("Proposed", 
                                "Multinomial MLE",
                                "Arbel et. al", 
                                "Chao & Shen", 
                                "iNEXT",
                                "Zero-replace"), 
                     values = c("Proposed" = "blue",
                                "Arbel et. al" = "#999999", 
                                "Chao & Shen" = "#E69F00",
                                "iNEXT" = "#56B4E9",
                                "Multinomial MLE" = "red",
                                "Zero-replace" = "black"))  +
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) 
bc_plot

euc_plot <- ggplot(full_data_frame, 
                   aes(x = q_char, 
                       y = mse_loss_euclidean, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab("Q: Number of taxa") + 
  ylab("MSE: Euclidean index") +
  theme_bw() + 
  scale_color_manual(breaks = c("Proposed", 
                                "Multinomial MLE",
                                "Arbel et. al", 
                                "Chao & Shen", 
                                "iNEXT",
                                "Zero-replace"), 
                     values = c("Proposed" = "blue",
                                "Arbel et. al" = "#999999", 
                                "Chao & Shen" = "#E69F00",
                                "iNEXT" = "#56B4E9",
                                "Multinomial MLE" = "red",
                                "Zero-replace" = "black"))  +
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) 
euc_plot

ggpubr::ggarrange(shannon_plot, 
                  simpson_plot, 
                  bc_plot, 
                  euc_plot, ncol=2, nrow=2, 
                  common.legend = T, legend="right")
ggsave("vary_q.pdf", 
       width=5.5, height = 3.5)
