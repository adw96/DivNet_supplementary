# Amy Willis, 26 November 2019
# A script to investigate the performance of divnet with LV models

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

source("../simulator_files/model_functions_lv.R")
source("../simulator_files/method_functions_lv.R")
source("../simulator_files/eval_functions.R")

#### Example of simulation

set.seed(191227)
eq1 <- exp(rnorm(20, 0, sd = 0.1))
long_run_abundances <- eq/sum(eq)
ta <- lv_draw(q = 20, ll = 60, 
              stoch_sd = 0.01,
              eq = eq1,  
              maxeig = 0.9)
ts1 <- ta$tc %>%
  as_tibble %>%
  bind_cols("time" = 1:nrow(ta$tc)) %>%
  mutate(time = time -1)%>%
  pivot_longer(cols = -time)
ts_plot <- ts1 %>%
  inner_join(tibble("name" = unique(ts1$name), 
                    "lr_value" = ta$eq)) %>%
  ggplot(aes(x = time, y = value, group = name, col = name,
             yintercept = lr_value)) +
  geom_line(lwd = 0.5) +
  theme_bw() + 
  ylab(expression('Absolute abundance V'[q]*'(t)'))+ 
  xlab("Time t") +
  geom_hline(aes(x = time, col = name,
                 yintercept = lr_value),
             lty = 2) + 
  xlim(0, nrow(ta$tc)) + 
  theme(legend.position = "none")
ts_plot

#### LV model, varying Q
lv_with_q_sim5 <- new_simulation(name = "lv-sim-1-longer-5",
                                 label = "LV with more taxa 5") %>%
  generate_model(make_lv, 
                 q = as.list(c(45, 30, 15)), 
                 ll = 80,
                 abundances_sd = 0.1,
                 stochasticity_sd = 0.01, 
                 mm = 1e5,
                 vary_along = c("q"), 
                 seed=191127) %>%
  simulate_from_model(nsim = 50) %>%
  run_method(list(divnet_method_lv,
                  plug_lv, 
                  inext_lv,
                  chao_shen_lv)) %>%
  evaluate(list(mse_loss_shannon_lv))

#### LV model, varying T
lv2_longer <- new_simulation(name = "lv-sim2-longer",
                             label = "LV with up to 100 time points") %>%
  generate_model(make_lv, 
                 q = 20, 
                 ll = as.list(c(20, 60, 100)),
                 abundances_sd = 0.1,
                 stochasticity_sd = 0.01, 
                 mm = 1e5,
                 vary_along = c("ll")) %>%
  simulate_from_model(nsim = 50)  %>%
  run_method(list(divnet_method_lv,
                  plug_lv, 
                  inext_lv,
                  chao_shen_lv)) %>%
  evaluate(list(mse_loss_shannon_lv))

###### Plots

lv_with_t_df <- right_join(lv2_longer %>% model %>% as.data.frame, 
                        lv2_longer %>% evals %>% as.data.frame, 
                        by = c("name" = "Model")) %>%
  as_tibble %>%
  mutate(Method = plyr::revalue(Method, 
                                c("divnet-lv"="Proposed",
                                  "chao-shen-lv" = "Chao &\nShen",
                                  "inext-lv" = "iNEXT",
                                  "plug-in-lv" = "Multinomial\nMLE"))) %>%
  mutate(t_char = paste(ifelse(ll == 100, "", " "), 
                        ll, " \n", Method , sep = "")) 
  # mutate(t_char = ifelse(substr(t_char) = 100, "", " "), 
  #                       ll, " \n", Method , sep = "")) 

lv_with_t <- lv_with_t_df %>%
  ggplot(aes(x = t_char, 
             y = mse_loss_shannon_lv, 
             col = Method)) +
  theme_bw() + scale_y_log10() + 
  geom_boxplot(outlier.size=0.1) + 
  xlab(label=expression("T: Number of timepoints")) +
  scale_color_manual(values=c("#56B4E9", "blue",  "black", "#999999"))  +
  ylab("MSE: Shannon diversity") + 
  theme(text = element_text(size = 12),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_discrete(breaks=lv_with_t_df$t_char,
                   labels=lv_with_t_df$t_char %>% substr(1,3))
lv_with_t

lv_with_q_df <- right_join(lv_with_q_sim5 %>% model %>% as.data.frame, 
                           lv_with_q_sim5 %>% evals %>% as.data.frame, 
                           by = c("name" = "Model")) %>%
  as_tibble %>%
  mutate(Method = plyr::revalue(Method, 
                                c("divnet-lv"="Proposed",
                                  "chao-shen-lv" = "Chao &\nShen",
                                  "inext-lv" = "iNEXT",
                                  "plug-in-lv" = "Multinomial\nMLE"))) %>%
  #mutate(q = q-1) %>%
  mutate(q_char = paste("", q, "\n", Method , sep = "")) 
lv_with_q <- lv_with_q_df %>%
  ggplot(aes(x = q_char, 
             y = mse_loss_shannon_lv, 
             col = Method)) +
  theme_bw() + scale_y_log10() + 
  geom_boxplot(outlier.size=0.1) + 
  xlab(label=expression("Q: Number of taxa")) +
  scale_color_manual(values=c("#56B4E9", "blue",  "black", "#999999"))  +
  ylab("MSE: Shannon diversity") + 
  theme(text = element_text(size = 12),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_x_discrete(breaks=lv_with_q_df$q_char,
                   labels=lv_with_q_df$q_char %>% substr(1,2))
lv_with_q

ggpubr::ggarrange(ts_plot, 
                  ggpubr::ggarrange(lv_with_q,
                                    lv_with_t,
                                    ncol=2, 
                                    common.legend = T, legend="bottom"), 
                  widths=c(0.75, 2))

ggpubr::ggarrange(ts_plot, 
                  ggpubr::ggarrange(lv_with_q,
                                    lv_with_t,
                                    ncol=2, 
                                    common.legend = T, legend="bottom"), 
                  nrow = 2, heights = c(0.33, 0.5))
ggsave("LV-robustness.pdf", units="in", width=8, height=6)

saveRDS(lv_with_q_df, "lv_with_q_df.RDS")
saveRDS(lv_with_t_df, "lv_with_t_df.RDS")
