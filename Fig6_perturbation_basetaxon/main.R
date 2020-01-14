# Amy Willis, 2 Jan 2020
# A script to create Figure 2 of the DivNet paper

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

## Set up data -- this is how Mike Lee constructed the data

# Begin Mike Lee's code

## this is an accompanying script of the workflow presented here: https://astrobiomike.github.io/amplicon/workflow_ex#analysis-in-r
## there is more detail presented if you follow along from the above page

### SETTING UP WORKING ENVIRONMENT; READING IN DATA
count_tab <- read.table("ASV_counts.txt", header=T, row.names=1, check.names=F)
tax_tab <- as.matrix(read.table("ASV_tax.txt", header=T, row.names=1, check.names=F, na.strings="", sep="\t"))
sample_info_tab <- read.table("sample_info.txt", header=T, row.names=1, check.names=F)

### TREATMENT OF "BLANKS"
# first we need to get a sum for each ASV across all 4 blanks and all 16 samples
blank_ASV_counts <- rowSums(count_tab[,1:4])
sample_ASV_counts <- rowSums(count_tab[,5:20])

# now we normalize them, here by dividing the samples' total by 4 – as there are 4x as many samples (16) as there are blanks (4)
norm_sample_ASV_counts <- sample_ASV_counts/4

# here we're getting which ASVs are deemed likely contaminants based on the threshold noted (http://localhost:4000/amplicon/workflow_ex#treatment-of-blanks):
blank_ASVs <- names(blank_ASV_counts[blank_ASV_counts * 10 > norm_sample_ASV_counts])
length(blank_ASVs) # this approach identified about 50 out of ~1550 that are likely to have orginated from contamination

# looking at the percentage of reads retained for each sample after removing these presumed contaminant ASVs shows that the blanks lost almost all of their sequences, while the samples, other than one of the bottom water samples, lost less than 1% of their sequences, as would be hoped
colSums(count_tab[!rownames(count_tab) %in% blank_ASVs, ]) / colSums(count_tab) * 100

# now that we've used our extraction blanks to identify ASVs that were likely due to contamination, we're going to trim down our count table by removing those sequences, and the blank samples, from further analysis
filt_count_tab <- count_tab[!rownames(count_tab) %in% blank_ASVs, -c(1:4)]
# make a filtered sample info table
filt_sample_info_tab<-sample_info_tab[-c(1:4), ]

# and let's add some colors to the sample info table that are specific to sample types and characteristics that we can use when plotting things
# we'll color the water samples blue: 
filt_sample_info_tab$color[filt_sample_info_tab$char == "water"] <- "blue"
# the biofilm sample a darkgreen:
filt_sample_info_tab$color[filt_sample_info_tab$char == "biofilm"] <- "darkgreen"
# the basalts with highly altered, thick outer rinds (>1 cm) brown ("chocolate4" is the best brown I can find...):
filt_sample_info_tab$color[filt_sample_info_tab$char == "altered"] <- "chocolate4"
# the basalts with smooth, glassy, thin exteriors black:
filt_sample_info_tab$color[filt_sample_info_tab$char == "glassy"] <- "black"
# and the calcified carbonate sample an ugly yellow:
filt_sample_info_tab$color[filt_sample_info_tab$char == "carbonate"] <- "darkkhaki"

filt_sample_info_tab

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

# model all the microbes in the rocks
n <- w_glassy_altered %>% sample_data %>% nrow 
my_w <- w_glassy_altered %>% otu_table %>% t %>% data.frame
my_w %>% dim
my_w %>% class
my_w %>% apply(1, sum)
which((my_w %>% apply(2, sum)) == 0)

my_x <- lm(rnorm(n) ~ w_glassy_altered %>% sample_data %$% char) %>% model.matrix
colnames(my_x) <- c("intercept", "glassy")
my_x

# ASV_2 is the most abundant across all samples, and has all non-zero entries
my_w %>% apply(2, function(x) all(x>0)) %>% which
my_w %>% apply(2, sum) %>% sort(decreasing=T) %>% head(5)
my_w[, "ASV_2"]
which(my_w %>% colnames == "ASV_2") # index 608 
tax_table(w_glassy_altered)["ASV_2",] # a Nitrospirae of Order Nitrospirales

1- mean(my_w > 0)
# dv_base <- readRDS(file="dv_base.RDS")
# dv_perturbations <- readRDS(file="dv_perturbations.RDS")
# perturbations <- c(1e-4, 1e-3, 1e-2, 0.05, seq(from = 0.1, to = 1, by = 0.1))


### PERTURBATIONS
set.seed(200102)
perturbations <- c(1e-4, 1e-3, 1e-2, 0.05, seq(from = 0.1, to = 1, by = 0.1))
dv_perturbations <- list()
for (i in 1:length(perturbations)) {
  dv_perturbations[[i]] <- DivNet::divnet(W = my_w,
                                          X = my_x,
                                          perturbation = perturbations[i],
                                          base = which(my_w %>% colnames == "ASV_2"), 
                                          ncores = 8,
                                          variance = "none")
}


####### BASE
set.seed(100) 
observed_in_all <- my_w %>% apply(2, function(x) all(x>0)) %>% which %>% names %>% sample(size=10) 
most_abundant <- my_w %>% apply(2, sum) %>% sort(decreasing=T) %>% head(10) %>% names
randomly_chosen <- names(my_w) %>% sample(size=10)
to_investigate <- c(observed_in_all, most_abundant, randomly_chosen) %>% unique %>% sort
to_investigate
dv_base <- list()
for (j in 1:length(to_investigate)) {
  dv_base[[j]] <- DivNet::divnet(W = my_w,
                                 X = my_x,
                                 perturbation = 0.5,
                                 base = which(my_w %>% colnames == to_investigate[j]), 
                                 ncores = 8,
                                 variance = "none")
}
dv_base


shannon_pert_plot <- tibble(perturbation = perturbations, 
                            Glassy = dv_perturbations %>%
                              lapply(function(x) x[[1]][[1]]$estimate) %>% unlist, # for glassy,
                            Altered = dv_perturbations %>%
                              lapply(function(x) x[[1]][[3]]$estimate) %>% unlist) %>%
  pivot_longer(cols=Glassy:Altered) %>%
  ggplot(aes(x = perturbation, y = value, col = name)) +
  geom_path() +
  ylim(0,10) + ylab("Estimate of Shannon index") + xlab("Perturbation parameter")+
  scale_colour_manual(values=c("#999999", "#56B4E9")) + 
  theme_bw() 
shannon_pert_plot

bc_pert_plot <- tibble(perturbation = perturbations, 
                       bc = dv_perturbations %>%
                         lapply(function(x) x[[3]][1,3]) %>% 
                         unlist) %>%
  ggplot(aes(x = perturbation, y = bc, col = "Glassy-Altered")) +
  geom_line() +
  ylim(0, 1) + 
  ylab("Estimate of Bray-Curtis index") + xlab("Perturbation parameter")+
  scale_colour_manual(values=c("purple")) + 
  theme_bw() 



shannon_base_plot <- tibble(base = to_investigate, 
                            Glassy = dv_base %>%
                              lapply(function(x) x[[1]][[1]]$estimate) %>% unlist, # for glassy,
                            Altered = dv_base %>%
                              lapply(function(x) x[[1]][[3]]$estimate) %>% unlist) %>%
  pivot_longer(cols=Glassy:Altered) %>%
  ggplot(aes(x = base, y = value, col = name)) +
  geom_point() +
  ylim(0,10) + ylab("Estimate of Shannon index") + xlab("Base taxon (unordered)")+
  scale_colour_manual(values=c("#999999", "#56B4E9")) + 
  theme_bw() 

bc_base_plot <- tibble(base = to_investigate, 
                       bc = dv_base %>%
                         lapply(function(x) x[[3]][1,3]) %>% 
                         unlist) %>%
  ggplot(aes(x = base, y = bc, col = "Glassy-Altered")) +
  geom_point() +
  ylim(0, 1) + 
  ylab("Estimate of Bray-Curtis index") + xlab("Base taxon (unordered)")+
  scale_colour_manual(values=c("purple")) + 
  theme_bw() 

### PLOTS

ggpubr::ggarrange(shannon_pert_plot + 
                    theme(legend.position="none") ,
                  shannon_base_plot + 
                    theme(legend.position="right", 
                                            legend.key.width=unit(2.7,"line"), 
                                            axis.text.x=element_blank()) + 
                    labs(col = ""), 
                  bc_pert_plot + 
                    theme(legend.position="none") ,
                  bc_base_plot + 
                    theme(legend.position="right", 
                                       legend.key.width=unit(0.1,"line"),
                                       axis.text.x=element_blank()) + 
                    labs(col = ""),
                  ncol=2, nrow=2, widths=c(1,1))
# ggsave("base_perturbation_parameter_glassy_altered.pdf", width=7, height = 5)

ggpubr::ggarrange(shannon_pert_plot + 
                    # theme(legend.position="none") +
                    NULL,
                  shannon_base_plot + 
                    theme(legend.position="right",
                          legend.key.width=unit(2.7,"line"),
                          axis.text.x=element_blank()) +
                    labs(col = ""), 
                  bc_pert_plot + 
                    # theme(legend.position="none") +
                    NULL,
                  bc_base_plot + 
                    theme(legend.position="right",
                          legend.key.width=unit(0.1,"line"),
                          axis.text.x=element_blank()) +
                    labs(col = ""),
                  ncol=4, nrow=1, common.legend=T)

shanns <- ggpubr::ggarrange(shannon_pert_plot + 
                    theme(legend.position="top",
                          legend.key.width=unit(2.7,"line"),
                          axis.text.x=element_blank(), 
                          axis.title.x=element_text(size=10))+
                      labs(col = ""),
                  shannon_base_plot + 
                    theme(legend.position="none",
                          axis.text.x=element_blank(), 
                          axis.title.x=element_text(size=10)),
                  ncol = 2, common.legend=T)

bcs <- ggpubr::ggarrange(bc_pert_plot + 
                    theme(legend.position="top",
                          legend.key.width=unit(2.7,"line"),
                          axis.text.x=element_blank(), 
                          axis.title.x=element_text(size=10))+
                      labs(col = ""),
                  bc_base_plot + 
                    theme(legend.position="none",
                          axis.text.x=element_blank(), 
                          axis.title.x=element_text(size=10)),
                  ncol=2, nrow=1, common.legend=T)
ggpubr::ggarrange(shanns, bcs, ncol=2)
# ggsave("base_perturbation_parameter_glassy_altered_long.pdf", width=8, height = 3)

## 
mean(my_w == 0) # percent zeroes
sum(apply(my_w > 0, 2, sum) == n)

## In Figure \ref{base_perturbation}, we observe sizeable changes in the diversity 
## estimates when varying $p$ close to zero (at most 26\%, -50\%, -24\% and -31\% 
## changes in Shannon, Simpson, Bray-Curtis, and Euclidean estimates for $\rho=0.001$ 
## compared to $\rho=0.5$), but smaller changes when $\rho$ is increased from 0.5 
## to 1 (at most 5\%, -24\%, -12\% and -13\% changes for $\rho=0.5$ to $\rho=1$).
sh_1e_3 <- dv_perturbations[perturbations == 0.001] %>%
  lapply(function(x) x[[1]][[1]]$estimate) %>% unlist
sh_0.5 <- dv_perturbations[perturbations == 0.5] %>%
  lapply(function(x) x[[1]][[1]]$estimate) %>% unlist # glassy
sh_1 <- dv_perturbations[perturbations == 1] %>%
  lapply(function(x) x[[1]][[1]]$estimate) %>% unlist # glassy

((sh_0.5-sh_1e_3)/sh_1e_3 * 100) %>% round # 26%
((sh_1-sh_0.5)/sh_0.5 * 100) %>% round # 5%

si_1e_3 <- dv_perturbations[perturbations == 0.001] %>%
  lapply(function(x) x[[2]][[1]]$estimate) %>% unlist
si_0.5 <- dv_perturbations[perturbations == 0.5] %>%
  lapply(function(x) x[[2]][[1]]$estimate) %>% unlist # glassy
si_1 <- dv_perturbations[perturbations == 1] %>%
  lapply(function(x) x[[2]][[1]]$estimate) %>% unlist # glassy
((si_0.5-si_1e_3)/si_1e_3 * 100) %>% round # -50%
((si_1-si_0.5)/si_0.5 * 100) %>% round # -24%

bc_1e_3 <- dv_perturbations[perturbations == 0.001] %>%
  lapply(function(x) x[[3]][1,3]) %>% unlist
bc_0.5 <- dv_perturbations[perturbations == 0.5] %>%
  lapply(function(x) x[[3]][1,3]) %>% unlist
bc_1 <- dv_perturbations[perturbations == 1] %>%
  lapply(function(x) x[[3]][1,3]) %>% unlist
((bc_0.5-bc_1e_3)/bc_1e_3 * 100) %>% round # -24%
((bc_1-bc_0.5)/bc_0.5 * 100) %>% round # -12%


eu_1e_3 <- dv_perturbations[perturbations == 0.001] %>%
  lapply(function(x) x[[4]][1,3]) %>% unlist
eu_0.5 <- dv_perturbations[perturbations == 0.5] %>%
  lapply(function(x) x[[4]][1,3]) %>% unlist
eu_1 <- dv_perturbations[perturbations == 1] %>%
  lapply(function(x) x[[4]][1,3]) %>% unlist
((eu_0.5-eu_1e_3)/eu_1e_3 * 100) %>% round # -31%
((eu_1-eu_0.5)/eu_0.5 * 100) %>% round # -13%

# saveRDS(dv_base, file="dv_base.RDS")
# saveRDS(dv_perturbations, file="dv_perturbations.RDS")
