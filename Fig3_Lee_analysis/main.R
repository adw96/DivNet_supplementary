# Amy Willis, 9 Jan 2020
# A script to create the data analysis figure

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

library(phyloseq)

## Same as for Fig 2:
## Set up data -- this is how Mike Lee constructed the data
## Begin Mike Lee's code -- see other script for detailed comments
count_tab <- read.table("../Fig6_perturbation_basetaxon/ASV_counts.txt", header=T, row.names=1, check.names=F)
tax_tab <- as.matrix(read.table("../Fig6_perturbation_basetaxon/ASV_tax.txt", header=T, row.names=1, check.names=F, na.strings="", sep="\t"))
sample_info_tab <- read.table("../Fig6_perturbation_basetaxon/sample_info.txt", header=T, row.names=1, check.names=F)
blank_ASV_counts <- rowSums(count_tab[,1:4])
sample_ASV_counts <- rowSums(count_tab[,5:20])
norm_sample_ASV_counts <- sample_ASV_counts/4
blank_ASVs <- names(blank_ASV_counts[blank_ASV_counts * 10 > norm_sample_ASV_counts])
length(blank_ASVs) # this approach identified about 50 out of ~1550 that are likely to have orginated from contamination
colSums(count_tab[!rownames(count_tab) %in% blank_ASVs, ]) / colSums(count_tab) * 100
filt_count_tab <- count_tab[!rownames(count_tab) %in% blank_ASVs, -c(1:4)]
filt_sample_info_tab<-sample_info_tab[-c(1:4), ]
filt_sample_info_tab$color[filt_sample_info_tab$char == "water"] <- "blue"
filt_sample_info_tab$color[filt_sample_info_tab$char == "biofilm"] <- "darkgreen"
filt_sample_info_tab$color[filt_sample_info_tab$char == "altered"] <- "chocolate4"
filt_sample_info_tab$color[filt_sample_info_tab$char == "glassy"] <- "black"
filt_sample_info_tab$color[filt_sample_info_tab$char == "carbonate"] <- "darkkhaki"
count_tab_phy <- otu_table(filt_count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)
sample_info_tab_phy <- sample_data(filt_sample_info_tab)
ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

## End Mike Lee's code
## Begin Amy Willis's code
# get only glassy and altered
w_glassy_altered <- ASV_physeq %>% 
  subset_samples(char %in% c("glassy", "altered")) %>%
  filter_taxa(function(x) sum(x >= 1) > 0, TRUE) # only taxa seen at least once in all remaining samples
n <- w_glassy_altered %>% sample_data %>% nrow 
my_w <- w_glassy_altered %>% otu_table %>% t %>% data.frame
my_x <- lm(rnorm(n) ~ w_glassy_altered %>% sample_data %$% char) %>% model.matrix
colnames(my_x) <- c("intercept", "glassy")
my_x


### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ###  Comparison ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### 

set.seed(200103)
fitted_model <- DivNet::fit_aitchison(W = my_w,
                                      X = my_x,
                                      perturbation = 0.5,
                                      base = which(my_w %>% colnames == "ASV_2"), 
                                      ncores = 4)
# saveRDS(fitted_model, "fitted_model.RDS")
sigma_eigenvalues <- fitted_model$sigma %>% eigen  
sigma_eigenvalues$values %>% sort %>% tail

plot(sigma_eigenvalues$values)
hist(sigma_eigenvalues$values)

12*20/60 # 4 hours


# 
# ### test from Johnstone, 2001 -- use library(RMTstat)
# l1 <- sigma_eigenvalues$values %>% sort %>% tail(1)
# p <- sigma_eigenvalues$values %>% length
# n <- nrow(my_w)
# # can reverse if p large -- see p300
# mu_np <- (sqrt(p-1) + sqrt(n))^2
# sigma_np <-  (sqrt(p-1) + sqrt(n))* (1/sqrt(p-1) + 1/sqrt(n))^(1/3)
# ts <- (l1-mu_np)/sigma_np
# c(l1, mu_np, sigma_np, ts)
# ts
# RMTstat::qtw(1e-5, beta = 1)
# x_dtw <- seq(from = -10, to = 10)
# dtw(x_dtw) 
# plot(x_dtw, dtw(x_dtw))


parametric_list <- list()
parametric_list_fitted <- list()
set.seed(1)
for (i in 1:20) {
  print(paste("i =", i))
  mw <- make_w(mu = fitted_model$fitted_y,
               Sigma = fitted_model$sigma,
               mm = apply(my_w, 1, sum), 
               base=which(my_w %>% colnames == "ASV_2"))
  parametric_list[[i]] <- DivNet::divnet(W = mw,
                                         X = my_x,
                                         perturbation = 0.5,
                                         base = which(my_w %>% colnames == "ASV_2"), 
                                         ncores = 8, 
                                         variance = "none")
}


# saveRDS(object=parametric_list, "parametric_list.RDS")
# parametric_list <- readRDS("parametric_list.RDS")


glassy_shannon_dv <- parametric_list %>% 
  lapply(function(x) x$shannon) %>%
  lapply(function(x) x[[1]])  %>% # for first obs; glassy
  lapply(function(x) x[["estimate"]]) %>% unlist
altered_shannon_dv <- parametric_list %>% 
  lapply(function(x) x$shannon) %>%
  lapply(function(x) x[[3]])  %>% # for first obs; altered
  lapply(function(x) x[["estimate"]]) %>% unlist
glassy_simpson_dv <- parametric_list %>% 
  lapply(function(x) x$simpson) %>%
  lapply(function(x) x[[1]])  %>% # for first obs; glassy
  lapply(function(x) x[["estimate"]]) %>% unlist
altered_simpson_dv <- parametric_list %>% 
  lapply(function(x) x$simpson) %>%
  lapply(function(x) x[[3]])  %>% # for first obs; altered
  lapply(function(x) x[["estimate"]]) %>% unlist


## Arbel
arbel_output <- DepGEM::gibbs(n.iter=500, Y = my_w, X = my_x[,2])

zz <- arbel_output$p_store 
get_shannon <- function(my_matrix) {
  my_matrix %>% apply(1, breakaway::true_shannon)
}
get_simpson <- function(my_matrix) {
  my_matrix %>% apply(1, breakaway::true_simpson)
}
arbel_bc <- apply(zz, 3, bray_curtis_true)
arbel_euclidean <- apply(zz, 3, euclidean_true)

apply(zz, 3, get_shannon) %>% c %>% plot # yep, converged

arbel_tibb <- tibble("Estimate" = apply(zz, 3, get_shannon)[my_x[,"glassy"] == 1, ] %>% c, 
                     "Sample" = "Glassy", 
                     "Estimator" = "Arbel",
                     "Estimand" = "Shannon") %>%
  bind_rows(tibble("Estimate" = apply(zz, 3, get_shannon)[my_x[,"glassy"] != 1, ] %>% c, 
                   "Sample" = "Altered", 
                   "Estimator" = "Arbel",
                   "Estimand" = "Shannon")) %>%
  bind_rows(tibble("Estimate" = apply(zz, 3, get_simpson)[my_x[,"glassy"] == 1, ] %>% c, 
                   "Sample" = "Glassy", 
                   "Estimator" = "Arbel",
                   "Estimand" = "Simpson")) %>%
  bind_rows(tibble("Estimate" = apply(zz, 3, get_simpson)[my_x[,"glassy"] != 1, ] %>% c, 
                   "Sample" = "Altered", 
                   "Estimator" = "Arbel",
                   "Estimand" = "Simpson")) %>%
  bind_rows(tibble("Estimate" = c(arbel_bc)[c(arbel_bc) > 0], 
                   "Sample" = "Glassy-Altered", 
                   "Estimator" = "Arbel",
                   "Estimand" = "Bray-Curtis")) %>%
  bind_rows(tibble("Estimate" = c(arbel_euclidean)[c(arbel_euclidean) > 0], 
                   "Sample" = "Glassy-Altered", 
                   "Estimator" = "Arbel",
                   "Estimand" = "Euclidean"))
# saveRDS(arbel_tibb, "arbel_tibb.RDS")

## Multinomial
df_x <- w_glassy_altered %>% 
  sample_data %>% 
  as.matrix %>% 
  as.data.frame %>% 
  rownames_to_column %>% 
  as_tibble
bc_w <- bray_curtis_true(my_w / rowSums(my_w))
euc_w <- euclidean_true(my_w / rowSums(my_w))
which_diff <- outer(my_x[ , "glassy"], 1-my_x[ , "glassy"])

mult_tibb <- w_glassy_altered %>% 
  sample_shannon %>% 
  summary %>%
  g %>%
  select(estimate, char) %>%
  rename(Estimate = estimate, Sample = char) %>%
  mutate("Estimator" = "Multinomial",
         "Estimand" = "Shannon") %>%
  bind_rows(w_glassy_altered %>% 
              sample_simpson %>% 
              summary %>%
              right_join(df_x, 
                         by = c("sample_names" = "rowname")) %>%
              select(estimate, char) %>%
              rename(Estimate = estimate, Sample = char) %>%
              mutate("Estimator" = "Multinomial",
                     "Estimand" = "Simpson")) %>%
  bind_rows(tibble("Estimate" = bc_w[which(which_diff==1)], 
                   "Sample" = "Glassy-Altered", 
                   "Estimator" = "Multinomial",
                   "Estimand" = "Bray-Curtis")) %>%
  bind_rows(tibble("Estimate" = euc_w[which(which_diff==1)], 
                   "Sample" = "Glassy-Altered", 
                   "Estimator" = "Multinomial",
                   "Estimand" = "Euclidean"))
mult_tibb



## iNEXT
reformat_draw <- my_w %>% t
rownames(reformat_draw) <- paste("Species", 1:(dim(reformat_draw)[1]), sep="")
colnames(reformat_draw) <- paste("Site", 1:(dim(reformat_draw)[2]), sep="")
my_table <- iNEXT::iNEXT(reformat_draw)$AsyEst
iNEXT_shannon <- data.frame("Estimate" = my_table[my_table[, "Diversity"] ==  "Shannon diversity", "Estimator"] %>% log, 
                            "Sample" = ifelse(my_x[,"glassy"] == 1, "Glassy", "Altered"),
                            "Estimator" = "iNEXT",
                            "Estimand" = "Shannon")
iNEXT_simpson <- data.frame("Estimate" = 1/my_table[my_table[, "Diversity"] == "Simpson diversity", "Estimator"], 
                            "Sample" = ifelse(my_x[,"glassy"] == 1, "Glassy", "Altered"),
                            "Estimator" = "iNEXT",
                            "Estimand" = "Simpson")



## Chao-Shen
f_tables <- apply(my_w, 1, breakaway::make_frequency_count_table)
f_tables2 <- f_tables %>% lapply(breakaway:::check_format)
ests <- f_tables2 %>% lapply(breakaway:::chao_shen_estimate) %>% unlist
chao_shen_shannon <- data.frame("Estimate" = ests, 
                                "Sample" = ifelse(my_x[,"glassy"] == 1, "Glassy", "Altered"),
                                "Estimator" = "Chao-Shen",
                                "Estimand" = "Shannon")


## Zero-replace
draw_imputed <- my_w
draw_imputed[which(draw_imputed == 0, arr.ind = TRUE)] <- 0.5
draw_imputed <- draw_imputed / apply(draw_imputed, 1, sum)

bc_w_zr <- bray_curtis_true(draw_imputed / rowSums(draw_imputed))
euc_w_zr <- euclidean_true(draw_imputed / rowSums(draw_imputed))

zr_tibb <- draw_imputed %>%
  apply(1, true_shannon) %>% 
  as_tibble %>%
  rename("Estimate" = value) %>%
  bind_cols(df_x) %>%
  select(Estimate, char) %>%
  rename(Sample = char) %>%
  mutate("Estimator" = "Zero-replace",
         "Estimand" = "Shannon") %>%
  bind_rows(draw_imputed %>% 
              apply(1, true_simpson) %>% 
              as_tibble %>%
              rename("Estimate" = value) %>%
              bind_cols(df_x) %>%
              select(Estimate, char) %>%
              rename(Sample = char) %>%
              mutate("Estimator" = "Zero-replace",
                     "Estimand" = "Simpson")) %>%
  bind_rows(tibble("Estimate" = bc_w_zr[which(which_diff==1)], 
                   "Sample" = "Glassy-Altered", 
                   "Estimator" = "Zero-replace",
                   "Estimand" = "Bray-Curtis")) %>%
  bind_rows(tibble("Estimate" = euc_w_zr[which(which_diff==1)], 
                   "Sample" = "Glassy-Altered", 
                   "Estimator" = "Zero-replace",
                   "Estimand" = "Euclidean"))



all_comparison_data <- arbel_tibb %>%
  bind_rows(mult_tibb) %>%
  bind_rows(iNEXT_shannon) %>%
  bind_rows(iNEXT_simpson) %>%
  bind_rows(chao_shen_shannon) %>% 
  bind_rows(zr_tibb) 
# write_csv(all_comparison_data, "all_comparison_data.csv")
# saveRDS(all_comparison_data, "all_comparison_data.RDS")

all_comparison_data_dv <- tibble("Estimate" = glassy_shannon_dv, 
                                 "Sample" = "glassy",
                                 "Estimator" = "Proposed",
                                 "Estimand" = "Shannon") %>%
  bind_rows(tibble("Estimate" = altered_shannon_dv, 
                   "Sample" = "altered",
                   "Estimator" = "Proposed",
                   "Estimand" = "Shannon")) %>%
  bind_rows(tibble("Estimate" = altered_simpson_dv, 
                   "Sample" = "altered",
                   "Estimator" = "Proposed",
                   "Estimand" = "Simpson")) %>%
  bind_rows(tibble("Estimate" = glassy_simpson_dv, 
                   "Sample" = "glassy",
                   "Estimator" = "Proposed",
                   "Estimand" = "Simpson")) %>%
  bind_rows(tibble("Estimate" = parametric_list %>% 
                     lapply(function(x) x$`bray-curtis`) %>%
                     lapply(function(x) x[1,3])  %>% unlist, 
                   "Sample" = "Glassy-Altered",
                   "Estimator" = "Proposed",
                   "Estimand" = "Bray-Curtis")) %>%
  bind_rows(tibble("Estimate" = parametric_list %>% 
                     lapply(function(x) x$`euclidean`) %>%
                     lapply(function(x) x[1,3])  %>% unlist, 
                   "Sample" = "Glassy-Altered",
                   "Estimator" = "Proposed",
                   "Estimand" = "Euclidean")) %>%
  bind_rows(all_comparison_data)
# saveRDS(all_comparison_data_dv, "all_comparison_data_dv.RDS")

all_comparison_data_dv %>% 
  mutate(Sample = tools::toTitleCase(Sample)) %>%
  mutate(Estimand = factor(Estimand, levels=c("Shannon", "Simpson", "Bray-Curtis", "Euclidean"))) %>% 
  ggplot(aes(x = Estimator, y = Estimate, col = Sample, fill = Sample)) +
  geom_boxplot(outlier.color="white", alpha = 0.5) +
  facet_wrap(~Estimand, scales="free") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("")  +
  scale_color_discrete(name = "Sample(s): ",
                       labels = c("Altered (n=8)", "Glassy (n=4)", "Glassy-Altered")) +
  scale_fill_discrete(name = "Sample(s): ",
                      labels = c("Altered (n=8)", "Glassy (n=4)", "Glassy-Altered")) +
  NULL

# ggsave("lee_dataset_analysis_glassy_altered_only.pdf", width=10, height = 7)
