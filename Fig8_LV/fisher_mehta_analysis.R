library(tidyverse)
library(magrittr)
# download.file("http://physics.bu.edu/~pankajm/Code/LIMITS-FINAL.nb", 
#               destfile="Fisher_Mehta.nb")
# notebook <- read_file("Fisher_Mehta.nb") 

# Steps
# - Take lines 7062 - 9868 from Fisher_Mehta.nb
# - Replace all `{` with `c(`
# - Replace all `}` with `)`
# - Replace all `"` with ``
# - Replace backticks with nothing
# - add a `my_data <- cbind(` at the start
# - save as fisher_mehta_lines_7062_9868.R

source("fisher_mehta_lines_7062_9868.R")
my_data %>% dim # 401 columns, 20 taxa
abundances <- my_data %>%  t
colnames(abundances) <- paste("taxon", 1:ncol(abundances), sep = "")
abundances %>%
  as_tibble %>%
  bind_cols("time" = 0:(-1 + nrow(abundances))) %>%
  pivot_longer(cols=taxon1:taxon20, names_to="taxon") %>%
  filter(taxon %in% paste("taxon", 1:6, sep = "")) %>%
  # filter(time %% )
  ggplot(aes(x = time, y = value, group = taxon, col = taxon)) +
  geom_line() +
  facet_wrap(~taxon)


relative_abundances <- abundances %>%
  as_tibble %>%
  bind_cols("time" = 0:(-1 + nrow(abundances)), .) %>%
  mutate_at(vars(matches("taxon")), function(x) {x/sum(x)})

relative_abundances
