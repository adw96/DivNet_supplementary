random_time_series <- function(q, ll, pp, stoch_sd, eq) {
  # simulates a timeseries of length ll for interaction matrix pp
  # and stochasticity stoch_sd
  
  # x: the stochasticity matrix (paper calls this eta)
  x <- matrix(data = exp(rnorm(q * (ll + 1), 
                               mean=0, 
                               sd=stoch_sd)),
              ncol = q, 
              nrow = ll + 1)
  tc <- matrix(NA, ncol=q, nrow = ll + 1)
  tc[1, ] <- x[1, ]
  for (tt in 1:ll) {
    tc[tt + 1, ] <- x[tt + 1, ] * tc[tt, ] * exp(pp %*% (tc[tt, ] - eq))
  }
  tc
}

lv_draw <- function(q, 
                    ll, 
                    stoch_sd, 
                    eq, 
                    maxeig = 0.9) {
  # rewriting the functions from http://physics.bu.edu/~pankajm/Code/LIMITS-FINAL.nb into R
  
  # diagonal interaction matrix
  p0 <- diag(runif(n = 1, min = -(1 + maxeig), max = (maxeig - 1)) / eq)
  
  # initialise
  possible_indices <- which(p0 == 0, arr.ind=TRUE)
  nonzeros <- (which(p0 != 0, arr.ind=TRUE) %>% nrow) - q
  it <- 1
  
  # randomly add in a new interaction
  while (nonzeros < 10 & it < (q*ll)^3) {
    # sod: selected_offdiagonals
    sod <- possible_indices[sample(1:nrow(possible_indices), size=1), ]
    
    # confirmed parameterisation of beta distn is the same between wolfram & R
    tmp <- (2*rbeta(n=1, shape1=1, shape2=1) - 1) * abs(diag(p0)[sod[1]]) 
    
    p1 <- p0
    p1[sod]
    p1[sod[1], sod[2]] <- tmp
    tc <- random_time_series(q = q, ll = ll, pp = p1, stoch_sd=stoch_sd, eq = eq) 
    if (!any(is.nan(tc))) {
      p0 <- p1
      possible_indices <- which(p0 == 0, arr.ind=TRUE)
      nonzeros <- (which(p0 != 0, arr.ind=TRUE) %>% nrow) - q
      
    } #else {
      # warning(print(tc))
    # }
    # else, try again
    it <- it + 1
  }
  list(pp = p0, tc = tc, its = it, eq = eq)
}

lv_draw_counts <- function(q, ll, stoch_sd, eq,  maxeig = 0.9, mm) {
  
  true_abundances <- lv_draw(q = q, ll = ll, 
                             stoch_sd = stoch_sd,
                             eq = eq,  maxeig = maxeig)
  
  relative_abundances_long <- true_abundances$tc %>%
    as_tibble %>%
    bind_cols("time" = 1:nrow(true_abundances$tc)) %>%
    pivot_longer(cols = -time) %>%
    group_by(time) %>%
    mutate(value = value/sum(value)) %>%
    ungroup 
  
  relative_abundances_wide <- relative_abundances_long %>%
    pivot_wider(names_from=name, values_from=value) %>%
    select(-time) %>% 
    as.matrix
  
  my_w <- foreach(i = 1:(ll+1), .combine = rbind) %do% {
    stats::rmultinom(n = 1, size = mm, prob = relative_abundances_wide[i,  ]) %>% c
  }
  
  my_w
}
make_lv <- function(q, 
                    ll,
                    abundances_sd = 0.1,
                    stochasticity_sd = 0.01, 
                    mm = 1e5) {
  
  # equilibrium abundances: `eq = RandomReal[LogNormalDistribution[0.0, 0.1], n];`
  eq <- exp(rnorm(q, 0, sd = abundances_sd))   
  long_run_abundances <- eq/sum(eq)
  
  new_model(name = "pLV",
            label = sprintf("Fisher & Mehta LV (q = %s, length = %s, sigma_abund = %s, sigma_stoch = %s)", 
                            q, ll, abundances_sd, stochasticity_sd),
            params = list(q = q, 
                          ll = ll,
                          zz = long_run_abundances, 
                          stochasticity_sd = stochasticity_sd,
                          abundances_sd = abundances_sd,
                          shannon = breakaway::true_shannon(long_run_abundances),
                          mm = mm),
            simulate = function(q, ll, zz, stochasticity_sd, mm, nsim) {
              # this function must return a list of length nsim
              replicate(nsim,
                        lv_draw_counts(q = q, 
                                       eq=zz,
                                       stoch_sd=stochasticity_sd, 
                                       ll = ll, 
                                       mm=mm), 
                        simplify=F)
              
            })
}

# library(foreach)
# lv_draw_counts(q = 5, 
#                stoch_sd=0.01, 
#                eq = exp(rnorm(5, 0, sd = 0.1)),
#                ll = 10, mm = 10000)
# 
# my_sim <- lv_draw(q = 5, 
#                   stoch_sd=0.01, 
#                   eq = exp(rnorm(5, 0, sd = 0.1)),
#                   ll = 10)
# my_sim$tc[,-1]/my_sim$tc[, 1]
# lv_draw_counts(q = 5, 
#                stoch_sd=0.01, 
#                eq = exp(rnorm(5, 0, sd = 0.1)),
#                ll = 10, mm = 100)
# my_sim$pp
# ts1 <- my_sim$tc %>%
#   as_tibble %>%
#   bind_cols("time" = 1:nrow(my_sim$tc)) %>% 
#   pivot_longer(cols = -time) 
# ts1 %>% 
#   inner_join(tibble("name" = unique(ts1$name), "lr_value" = my_sim$eq)) %>%
#   # group_by(time) %>%
#   # mutate(value = value/sum(value)) %>%
#   ggplot(aes(x = time, y = value, group = name, col = name, 
#              yintercept = lr_value)) +
#   geom_line() +
#   geom_hline(aes(x = time, col = name, 
#                  yintercept = lr_value), 
#              lty = 2)

