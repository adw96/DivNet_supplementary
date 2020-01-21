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



lv1 <- new_simulation(name = "lv-sim-1",
                      label = "How does DivNet do with LV models?") %>%
  generate_model(make_lv, 
                 q = as.list(c(5, 10, 15, 20)), 
                 ll = 20,
                 abundances_sd = 0.1,
                 stochasticity_sd = 0.01, 
                 mm = 1e5,
                 vary_along = c("q")) %>%
  simulate_from_model(nsim = 50) %>%
  run_method(list(divnet_method_lv,
                  plug_lv, 
                  inext_lv,
                  chao_shen_lv)) %>%
  evaluate(list(mse_loss_shannon_lv))

plot_eval_by(lv1, metric_name="mse_loss_shannon_lv", varying="q")

right_join(lv1 %>% model %>% as.data.frame, 
           lv1 %>% evals %>% as.data.frame, 
           by = c("name" = "Model")) %>%
  ggplot(aes(x = interaction(Method, q), 
             y = mse_loss_shannon_lv, 
             col = Method)) +
  geom_boxplot() 
# median better for divnet for all q
# spread potentially worse with q


lv2 <- new_simulation(name = "lv-sim2",
                      label = "LV with l") %>%
  generate_model(make_lv, 
                 q = 20, 
                 ll = as.list(c(10, 20, 30, 40)),
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

right_join(lv2 %>% model %>% as.data.frame, 
           lv2 %>% evals %>% as.data.frame, 
           by = c("name" = "Model")) %>%
  ggplot(aes(x = interaction(Method, ll), 
             y = mse_loss_shannon_lv, 
             col = Method)) +
  geom_boxplot() + 
  scale_y_sqrt()
# median better for divnet for all q
# spread better with ll


lv3 <- new_simulation(name = "lv-sim3",
                      label = "LV with abundances_sd") %>%
  generate_model(make_lv, 
                 q = 20, 
                 ll = 20,
                 abundances_sd = as.list(seq(0.05, 0.35, length.out = 4)),
                 stochasticity_sd = 0.01, 
                 mm = 1e5,
                 vary_along = c("abundances_sd")) %>%
  simulate_from_model(nsim = 50)  %>%
  run_method(list(divnet_method_lv,
                  plug_lv, 
                  inext_lv,
                  chao_shen_lv)) %>%
  evaluate(list(mse_loss_shannon_lv))
right_join(lv3 %>% model %>% as.data.frame, 
           lv3 %>% evals %>% as.data.frame, 
           by = c("name" = "Model")) %>%
  ggplot(aes(x = interaction(Method, abundances_sd), 
             y = mse_loss_shannon_lv, 
             col = Method)) +
  geom_boxplot()
right_join(lv3 %>% model %>% as.data.frame, 
           lv3 %>% evals %>% as.data.frame, 
           by = c("name" = "Model")) %>%
  ggplot(aes(x = interaction(Method, abundances_sd), 
             y = mse_loss_shannon_lv, 
             col = Method)) +
  geom_boxplot() 
# # better with higher abundances_sd

# right_join(lv3 %>% model %>% as.data.frame, 
#            lv3 %>% evals %>% as.data.frame, 
#            by = c("name" = "Model")) %>%
#   mutate(shannon = round(shannon, 3)) %>% 
#   ggplot(aes(x = interaction(Method, shannon), 
#              y = mse_loss_shannon_lv, 
#              col = Method)) +
#   geom_boxplot()

### test everything works

test_simulation <- new_simulation(name = "lv-sim-test",
                                  label = "LV test") %>%
  generate_model(make_lv, 
                 q = as.list(c(5, 10)), 
                 ll = 20,
                 abundances_sd = 0.1,
                 stochasticity_sd = 0.01, 
                 mm = 1e5,
                 vary_along = c("q")) %>%
  simulate_from_model(nsim = 3) 
source("../simulator_files/method_functions_lv.R")
test_simulation %>%
  run_method(list(divnet_method_lv,
                  plug_lv, 
                  inext_lv,
                  chao_shen_lv))
m <- model(test_simulation, subset = "pLV/abundances_sd_0.1/ll_20/mm_100000/q_5/stochasticity_sd_0.01")
d <- draws(test_simulation, subset = "pLV/abundances_sd_0.1/ll_20/mm_100000/q_5/stochasticity_sd_0.01", index = 1)
chao_shen_lv@method(model = m, draw = d@draws$r1.1)


### 


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


lv1_longer <- new_simulation(name = "lv-sim-1-longer",
                             label = "LV with more taxa") %>%
  generate_model(make_lv, 
                 q = as.list(c(20, 40, 60, 80, 100)), 
                 ll = 60,
                 abundances_sd = 0.1,
                 stochasticity_sd = 0.01, 
                 mm = 1e5,
                 vary_along = c("q")) %>%
  simulate_from_model(nsim = 50) %>%
  run_method(list(divnet_method_lv,
                  plug_lv, 
                  inext_lv,
                  chao_shen_lv)) %>%
  evaluate(list(mse_loss_shannon_lv))

lv_with_q <- right_join(lv1_longer %>% model %>% as.data.frame, 
                        lv1_longer %>% evals %>% as.data.frame, 
                        by = c("name" = "Model")) %>%
  as_tibble %>%
  mutate(Method = plyr::revalue(Method, 
                                c("divnet-lv"="Proposed",
                                  "chao-shen-lv" = "Chao & Shen",
                                  "inext-lv" = "iNEXT",
                                  "plug-in-lv" = "Multinomial MLE"))) %>%
  filter(q %in% c(20, 60, 80)) %>%
  # filter(q %in% c(20, 40, 60, 80)) %>%
  mutate(q_char = paste("Q = ", q, "\n", Method , sep = "")) %>%
  ggplot(aes(x = q_char, 
             y = mse_loss_shannon_lv, 
             col = Method)) +
  theme_bw() + 
  scale_y_log10() + 
  geom_boxplot(outlier.size=0.1) + 
  xlab(label=expression("Q: Number of taxa")) +
  scale_color_manual(values=c("#56B4E9", "blue",  "black", "#999999"))  +
  ylab("MSE: Shannon diversity") + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
lv_with_q 


lv_with_t <- right_join(lv2_longer %>% model %>% as.data.frame, 
                        lv2_longer %>% evals %>% as.data.frame, 
                        by = c("name" = "Model")) %>%
  as_tibble %>%
  mutate(Method = plyr::revalue(Method, 
                                c("divnet-lv"="Proposed",
                                  "chao-shen-lv" = "Chao & Shen",
                                  "inext-lv" = "iNEXT",
                                  "plug-in-lv" = "Multinomial MLE"))) %>%
  # filter(ll %in% c(20, 80)) %>%
  filter(ll %in% c(20, 40, 60, 80, 100)) %>%
  mutate(ll = ll-1) %>%
  mutate(q_char = paste("T = ", ll, "\n", Method , sep = "")) %>%
  ggplot(aes(x = q_char, 
             y = mse_loss_shannon_lv, 
             col = Method)) +
  theme_bw() + 
  scale_y_log10() + 
  geom_boxplot(outlier.size=0.1) + 
  xlab(label=expression("T: Number of timepoints")) +
  scale_color_manual(values=c("#56B4E9", "blue",  "black", "#999999"))  +
  ylab("MSE: Shannon diversity") + 
  theme(text = element_text(size = 9),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
lv_with_t

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
  pivot_longer(cols = -time)
ts_plot <- ts1 %>%
  inner_join(tibble("name" = unique(ts1$name), "lr_value" = ta$eq)) %>%
  ggplot(aes(x = time, y = value, group = name, col = name,
             yintercept = lr_value)) +
  geom_line(lwd = 0.5) +
  theme_bw() + 
  geom_hline(aes(x = time, col = name,
                 yintercept = lr_value),
             lty = 2) 
ts_plot

ggpubr::ggarrange(ts_plot, 
                  ggpubr::ggarrange(lv_with_q,
                                    lv_with_t,
                                    ncol=2, 
                                    common.legend = T, legend="right"), widths=c(1,2))
