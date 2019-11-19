# Amy Willis, 18 November 2019
# A script to create Figure 2 of the DivNet paper
# n vs MSE

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
source("../simulator_files/eval_functions.R")

## @knitr main
sim_alpha_sigma <- new_simulation(name = "alpha-diversity-estimation-sigma",
                                  label = "Estimating alpha diversity under Aitchison models over varying sigma_max") %>%
  generate_model(make_aitchison, seed = 1234,
                 n = 20, #as.list(c(10, 28, 46, 62, 80)),
                 p = 2, 
                 q = 20,
                 sigma_min = 0.01,
                 sigma_max = as.list(seq(from = 0.01, to = 10, length.out = 5)),
                 beta_sd = 1,
                 vary_along = c("sigma_max")) %>%
  simulate_from_model(nsim = 5) #100)

source("method_functions.R")
sim_alpha_sigma_methods <- sim_alpha_sigma %>%
  run_method(list(chao_shen,
                  divnet_method,
                  arbel,
                  inext,
                  plug), 
             parallel = list(socket_names = 6))

sim_alpha_sigma_evaluated <- sim_alpha_sigma_methods %>% ## if want in parallel
  evaluate(list(mse_loss_shannon, mse_loss_simpson, 
                mse_loss_euclidean, mse_loss_bray_curtis, 
                mse_loss_bray_curtis_subset, mse_loss_euclidean_subset)) 

# sim_alpha <- load_simulation("alpha-diversity-estimation") 

## @knitr plots
# save_simulation(sim_alpha_sigma_evaluated)

simpson_plot <- plot_eval_by(sim_alpha_sigma_evaluated, "mse_loss_simpson", varying = "sigma_max")
shannon_plot <- plot_eval_by(sim_alpha_sigma_evaluated, "mse_loss_shannon", varying = "sigma_max")
time_plot <- plot_eval_by(sim_alpha_sigma_evaluated, "time", varying = "sigma_max")

#pdf("vary_sigma_max_beta_sd_1_n_20_p_2_q_20_sigma_min_0_01.pdf")
grid.arrange(simpson_plot, shannon_plot, time_plot, nrow = 3)
dev.off()
## @knitr tables

tabulate_eval(sim, "mse_loss_shannon", #output_type = "markdown",
              format_args = list(digits = 1))


###############################################################
########################## alpha diversity -- n ###############
###############################################################

source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")
sim_alpha_n <- new_simulation(name = "alpha-diversity-estimation-with-n",
                              label = "Estimating alpha diversity over varying n") %>%
  generate_model(make_aitchison, seed = 12345,
                 n = as.list(c(40, 30, 20, 10)),
                 p = 2, 
                 q = 20,
                 sigma_min = 0.01,
                 sigma_max = 5,
                 beta_sd = 1,
                 vary_along = c("n")) %>%
  simulate_from_model(nsim = 100)  

sim_alpha_n_methods <- sim_alpha_n %>%
  run_method(list(chao_shen,
                  #aitchison_poor, 
                  arbel,
                  inext,
                  plug), parallel = list(socket_names = 6)) 

#save_simulation(sim_alpha_n_methods)
sim_alpha_n_methods %>% ## if want in parallel
  evaluate(list(mse_loss_shannon, mse_loss_simpson, 
                mse_loss_euclidean, mse_loss_bray_curtis, 
                mse_loss_bray_curtis_subset, mse_loss_euclidean_subset)) -> sim_alpha_n


sim_alpha_n_methods2 <- sim_alpha_n_methods  %>% 
  simulator::rename("n-add-more") %>%
  relabel("Effect of varying n") %>%
  run_method(list(plug,
                  plug_add_half,
                  plug_pooled), 
             parallel = list(socket_names = 6))  %>% 
  evaluate(list(mse_loss_shannon, mse_loss_simpson, 
                mse_loss_euclidean, mse_loss_bray_curtis, 
                mse_loss_bray_curtis_subset, mse_loss_euclidean_subset))

evals(sim_alpha_n_methods2) %>% as.data.frame

# my_data_frame <- evals(sim_alpha_n) %>% as.data.frame
my_data_frame <- evals(sim_alpha_n_methods2) %>% as.data.frame

ns <- my_data_frame$Model %>% as.character %>% strsplit("/") %>% lapply(function(x) x[3]) %>% 
  unlist %>% strsplit(split="_") %>% lapply(function(x) x[2]) %>% unlist %>% as.numeric

full_data_frame <- data.frame(my_data_frame, "n" = ns, "n_char" = as.factor(ns),
                              "n_method" = paste(ns, "_", my_data_frame$Method, sep = ""))
full_data_frame$Method %>% unique
full_data_frame$Method <- plyr::revalue(full_data_frame$Method, c("aitchison-poor"="Proposed",
                                                                  "arbel" = "Arbel et. al",
                                                                  "chao-shen" = "Chao & Shen",
                                                                  "inext" = "iNEXT",
                                                                  "plug-in" = "Multinomial MLE",
                                                                  "plug-in-add-half" = "Zero-replace",
                                                                  "plug-in-pooled" = "Pooled Multinomial"))
full_data_frame$Method %<>% relevel("Multinomial MLE")
full_data_frame$Method %<>% relevel("Zero-replace")
full_data_frame$Method %<>% relevel("Pooled Multinomial")
full_data_frame$Method %<>% relevel("Proposed")

full_data_frame %<>% filter(Method != "Pooled Multinomial")

shannon_plot <- ggplot(full_data_frame[order(full_data_frame$n), ], 
                       aes(x = n_char, y = mse_loss_shannon, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab("n: Number of samples") +
  ylab("MSE: Shannon index") +
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00", "#56B4E9",  "black"))  
shannon_plot

simpson_plot <- ggplot(full_data_frame[order(full_data_frame$n), ], 
                       aes(x = n_char, y = mse_loss_simpson, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab("n: Number of samples") +
  ylab("MSE: Simpson index")+
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00", "black"))  
simpson_plot

bc_plot <- ggplot(full_data_frame[order(full_data_frame$n), ], 
                  aes(x = n_char, y = mse_loss_bray_curtis, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab("n: Number of samples") +
  ylab("MSE: Bray-Curtis index") +
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00", "#56B4E9"))  
bc_plot

euc_plot <- ggplot(full_data_frame[order(full_data_frame$n), ], 
                   aes(x = n_char, y = mse_loss_euclidean_subset, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab("n: Number of samples") +
  ylab("MSE: Euclidean index")+
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00", "#56B4E9")) 
euc_plot

time_plot <- ggplot(full_data_frame[order(full_data_frame$n), ], 
                    aes(x = n_char, y = time, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab("n: Number of samples") + 
  ylab("Computation time\n(seconds)")+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00", "#56B4E9",  "black"))  
time_plot

time_plot
ggsave("vary_n_beta_sd_1_sigma_max_5_p_2_q_20_sigma_min_0_01_time.pdf", height=2, width = 5.5)

ggpubr::ggarrange(shannon_plot, simpson_plot, bc_plot, euc_plot, ncol=2, nrow=2, common.legend = T, legend="right")
ggsave("vary_n_beta_sd_1_sigma_max_5_p_2_q_20_sigma_min_0_01_all_shared.pdf", width=5.5, height = 3.5)




###########################################################################
########################## alpha diversity -- sigma_max ###############
###########################################################################
source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")
sim_alpha_sigma_max <- new_simulation(name = "alpha-diversity-estimation-with-sigma-max",
                                      label = "Estimating alpha diversity over varying sigma_max") %>%
  generate_model(make_aitchison, seed = 123456,
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
                  aitchison_poor, 
                  arbel,
                  inext,
                  plug), parallel = list(socket_names = 6)) 

sim_alpha_sigma_max_methods2 <- sim_alpha_sigma_max_methods
sim_alpha_sigma_max_methods2 <- sim_alpha_sigma_max %>%
  run_method(list(plug), parallel = list(socket_names = 6)) 


#save_simulation(sim_alpha_sigma_max_methods)
sim_alpha_sigma_max_methods2 %>% 
  evaluate(list(mse_loss_shannon, mse_loss_simpson, 
                mse_loss_euclidean, mse_loss_bray_curtis, 
                mse_loss_bray_curtis_subset, mse_loss_euclidean_subset)) -> sim_alpha_sigma_max_done2

sim_alpha_sigma_max_methods3 <- sim_alpha_sigma_max_methods  %>% 
  simulator::rename("sigma-max-add-more") %>%
  relabel("Effect of varying sigma_max") %>%
  run_method(list(plug_add_half), 
             parallel = list(socket_names = 6))  %>% 
  evaluate(list(mse_loss_shannon, mse_loss_simpson, 
                mse_loss_euclidean, mse_loss_bray_curtis, 
                mse_loss_bray_curtis_subset, mse_loss_euclidean_subset))


# save_simulation(sim_alpha_sigma_max_done)
my_data_frame <- evals(sim_alpha_sigma_max_done) %>% as.data.frame
my_data_frame <- evals(sim_alpha_sigma_max_methods3) %>% as.data.frame

sms <- my_data_frame$Model %>% as.character %>% strsplit("/") %>% lapply(function(x) x[6]) %>% 
  unlist %>% strsplit(split="_") %>% lapply(function(x) x[3]) %>% unlist %>% as.numeric

full_data_frame <- data.frame(my_data_frame, "sigma_max" = sms, "sigma_max_char" = as.factor(sms),
                              "sigma_max_method" = paste(sms, "_", my_data_frame$Method, sep = ""))
full_data_frame$Method %>% unique
full_data_frame %<>% filter(Method != "plug-in-pooled")
full_data_frame$Method <- plyr::revalue(full_data_frame$Method, c("aitchison-poor"="Proposed",
                                                                  "arbel" = "Arbel et. al",
                                                                  "chao-shen" = "Chao & Shen",
                                                                  "inext" = "iNEXT",
                                                                  "plug-in" = "Multinomial MLE",
                                                                  "plug-in-add-half" = "Zero-replace"))

full_data_frame$Method %<>% relevel("Multinomial MLE")
full_data_frame$Method %<>% relevel("Zero-replace")
full_data_frame$Method %<>% relevel("Proposed")

shannon_plot <- ggplot(full_data_frame[order(full_data_frame$sigma_max), ], 
                       aes(x = sigma_max_char, y = mse_loss_shannon, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab(label=expression(sigma[max]~": Strength of cooccurrences")) +
  ylab("MSE: Shannon index")  +
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00", "#56B4E9",  "black"))  
shannon_plot

simpson_plot <- ggplot(full_data_frame[order(full_data_frame$sigma_max), ], 
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

bc_plot <- ggplot(full_data_frame[order(full_data_frame$sigma_max), ], 
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

euc_plot <- ggplot(full_data_frame[order(full_data_frame$sigma_max), ], 
                   aes(x = sigma_max_char, y = mse_loss_euclidean_subset, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab(label=expression(sigma[max]~": Strength of cooccurrences")) +
  ylab("Log MSE: Euclidean index") +
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00", "#56B4E9")) 
euc_plot

time_plot <- ggplot(full_data_frame[order(full_data_frame$sigma_max), ], 
                    aes(x = sigma_max_char, y = time, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab(label=expression(sigma[max]~": Strength of cooccurrences")) +
  ylab("Mean computational time (seconds)")+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00", "#56B4E9", "black")) 
time_plot

ggpubr::ggarrange(shannon_plot, simpson_plot, bc_plot, euc_plot, ncol=2, nrow=2, 
                  common.legend = T, legend="right")

# ggsave("vary_sigma_max_beta_sd_1_n_20_p_2_q_20_sigma_min_0_01_all_shared.pdf", width=6, height = 4)
# save.image("end180309.RData")
###########################################################################
########################## alpha diversity -- q ###############
###########################################################################
source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")
sim_q <- new_simulation(name = "diversity-estimation-with-q",
                        label = "Estimating diversity over varying q") %>%
  generate_model(make_aitchison, seed = 1234567,
                 n = 20,
                 p = 2, 
                 q = as.list(c(20, 60, 100)),
                 sigma_min = 0.01,
                 sigma_max = 5,
                 beta_sd = 1,
                 vary_along = c("q")) %>%
  simulate_from_model(nsim = 100)  

sim_q_methods <- sim_q %>%
  run_method(list(chao_shen,
                  aitchison_poor, 
                  arbel,
                  inext,
                  plug), parallel = list(socket_names = 6)) 
sim_q_methods %>% 
  evaluate(list(mse_loss_shannon, mse_loss_simpson, 
                mse_loss_euclidean, mse_loss_bray_curtis, 
                mse_loss_bray_curtis_subset, mse_loss_euclidean_subset)) #-> sim_q_done

## whoops!
sim_q_methods <- sim_q %>%
  run_method(list(plug_pooled, plug_add_half), parallel = list(socket_names = 6)) 


sim_q_methods %>% evaluate(list(mse_loss_shannon, mse_loss_simpson, 
                                mse_loss_euclidean, mse_loss_bray_curtis, 
                                mse_loss_bray_curtis_subset, mse_loss_euclidean_subset)) -> sim_q_methods_new_done

plot_eval(sim_q_methods_new_done, "mse_loss_shannon")
evals(sim_q_methods_new_done) %>% as.data.frame

#save_simulation(sim_alpha_sigma_max_methods)


# save_simulation(sim_alpha_sigma_max_done)
my_data_frame <- rbind(evals(sim_q_done) %>% as.data.frame,
                       evals(sim_q_methods_new_done) %>% as.data.frame)

sms <- my_data_frame$Model %>% as.character %>% strsplit("/") %>% lapply(function(x) x[5]) %>% 
  unlist %>% strsplit(split="_") %>% lapply(function(x) x[2]) %>% unlist %>% as.numeric

full_data_frame <- data.frame(my_data_frame, "q" = sms, "q_char" = as.factor(sms),
                              "q_method" = paste(sms, "_", my_data_frame$Method, sep = ""))

full_data_frame$Method %>% unique
full_data_frame %<>% filter(Method != "plug-in-pooled")
full_data_frame$Method %>% unique

full_data_frame$Method <- plyr::revalue(full_data_frame$Method, c("aitchison-poor"="Proposed",
                                                                  "arbel" = "Arbel et. al",
                                                                  "plug-in" = "Multinomial MLE",
                                                                  "plug-in-add-half" = "Zero-replace",
                                                                  "chao-shen" = "Chao & Shen",
                                                                  "inext" = "iNEXT"))

full_data_frame$Method %<>% relevel("Multinomial MLE")
full_data_frame$Method %<>% relevel("Zero-replace")
full_data_frame$Method %<>% relevel("Proposed")
shannon_plot <- ggplot(full_data_frame[order(full_data_frame$q), ], 
                       aes(x = q_char, y = mse_loss_shannon, col = Method)) +
  geom_boxplot(outlier.size = 0.1) + 
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Q: Number of taxa") +
  ylab("MSE: Shannon index")+
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00", "#56B4E9", "black")) 
shannon_plot

simpson_plot <- ggplot(full_data_frame[order(full_data_frame$q), ], 
                       aes(x = q_char, y = mse_loss_simpson, col = Method)) +
  geom_boxplot(outlier.size = 0.1) + 
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Q: Number of taxa") +
  ylab("MSE: Simpson index")  +
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00", "black")) 
simpson_plot

bc_plot <- ggplot(full_data_frame[order(full_data_frame$q), ], 
                  aes(x = q_char, y = mse_loss_bray_curtis, col = Method)) +
  geom_boxplot(outlier.size = 0.1) + 
  xlab("Q: Number of taxa") +
  ylab("MSE: Bray-Curtis index")+
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00")) 
bc_plot

euc_plot <- ggplot(full_data_frame[order(full_data_frame$q), ], 
                   aes(x = q_char, y = mse_loss_euclidean_subset, col = Method)) +
  geom_boxplot(outlier.size = 0.1) + 
  xlab("Q: Number of taxa") +
  ylab("MSE: Euclidean index") +
  theme_bw() + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00")) 
euc_plot

time_plot <- ggplot(full_data_frame[order(full_data_frame$q), ], 
                    aes(x = q_char, y = time, col = Method)) +
  geom_boxplot() + 
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #scale_y_log10() +
  xlab("Q: Number of taxa") +
  ylab("Mean computational time (seconds)") + 
  scale_color_manual(values=c("blue", "red", "#999999", "#E69F00", "#56B4E9",  "black", "purple"))
time_plot

ggpubr::ggarrange(shannon_plot, simpson_plot, bc_plot, euc_plot, ncol=2, nrow=2, common.legend = T, legend="right")

# ggsave("vary_q_beta_sd_1_n_20_p_2_q_20_sigma_min_0_01_sigma_max_5_all_shared.pdf", width=6, height = 3)

# save.image("end170409.R")
