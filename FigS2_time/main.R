# 

library(tidyverse)
n_df <- read_csv("../Fig2_MSE_n/n_df.csv")
q_df <- read_csv("../Fig5_MSE_numtaxa/q_df.csv")

n_df$Method %>% unique
n_df2 <- n_df %>%
  mutate(Method = factor(Method, levels = c("Proposed", 
                                            "Multinomial MLE",
                                            "Arbel et. al", 
                                            "Chao & Shen", 
                                            "iNEXT",
                                            "Zero-replace"))) %>%
  # mutate(Method = relevel(Method, "Multinomial MLE")) %>%
  # mutate(Method = relevel(Method, "Zero-replace")) %>%
  # mutate(Method = relevel(Method, "Proposed")) %>%
  # mutate(Method = plyr::revalue(Method, 
  #                               c("arbel" = "Arbel et. al",
  #                                 "divnet"="Proposed",
  #                                 "chao-shen" = "Chao &\nShen",
  #                                 "inext" = "iNEXT",
  #                                 "plug-in" = "Multinomial\nMLE",
  #                                 "plug-in-add-half" = "Zero-replace"))) %>%
  mutate(n_method = paste(n_char, " \n", Method , sep = "")) 

q_df2 <- q_df %>%
  mutate(q_method = paste(ifelse(q_char == 100, "", " "),
                          q_char, " \n", Method , sep = "")) %>%
  mutate(Method = as.factor(Method)) %>%
  mutate(Method = relevel(Method, "Multinomial MLE")) %>%
  mutate(Method = relevel(Method, "Zero-replace")) %>%
  mutate(Method = relevel(Method, "Proposed")) %>%
  mutate(q_char = factor(q_char, levels=c("20", "60", "100"))) 

time_plot_n <- n_df2 %>%
  filter(n_char > 10) %>%
  ggplot(aes(x = n_method, 
             y = time, col = Method)) +
  geom_boxplot(outlier.size=0.1) + 
  xlab("n: Number of samples") + 
  ylab("Computation time\n(seconds)") + 
  ylim(0, 80) +
  theme_bw() + 
  scale_x_discrete(breaks=n_df2$n_method,
                   labels=n_df2$n_method %>% substr(1,3)) +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 4) ,
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
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
                                "Zero-replace" = "black"))  
time_plot_n

time_plot_q <- ggplot(q_df2, 
                      aes(x = q_method, 
                          y = time, col = Method)) +
  geom_boxplot(outlier.size = 0.1) + 
  xlab("Q: Number of taxa") +
  ylab("Computation time\n(seconds)") + 
  theme_bw() + 
  ylim(0, 80) +
  scale_x_discrete(breaks=q_df2$q_method,
                   labels=q_df2$q_method %>% substr(1,3)) +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 4) ,
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
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
  NULL
time_plot_q

ggpubr::ggarrange(time_plot_n, time_plot_q, 
                  common.legend=TRUE, legend="right")
ggsave("timeplot.pdf", units="in", width=7, height=2)

