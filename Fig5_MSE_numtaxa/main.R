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
