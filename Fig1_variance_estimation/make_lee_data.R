# This is the main simulator file
source("/Users/adwillis/research/alpha-aitchison/data/lee/mike_phyloseq_object.R")
count_tab_phy <- otu_table(filt_count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)
sample_info_tab_phy <- sample_data(filt_sample_info_tab)

ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)
mike <- ASV_physeq
# otu_table()   OTU Table:         [ 1490 taxa and 16 samples ]
# sample_data() Sample Data:       [ 16 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 1490 taxa by 7 taxonomic ranks ]


## We only have 14 samples, so need fewer taxa
duplicates <- mike %>% subset_samples(char %in% c("glassy", "carbonate", "altered"))

duplicate_indices <- mike %>% sample_data %$% char %in% c("glassy", "carbonate", "altered")

# model all the microbes in the rocks
n <- duplicates %>% sample_data %$% char %>% length
my_w <- duplicates %>% otu_table %>% t %>% data.frame
my_w %>% dim
my_w %>% class

my_x <- lm(rnorm(n) ~ duplicates %>% sample_data %$% char) %>% model.matrix
colnames(my_x) <- c("intercept", "carbonate", "glassy")
my_x
