## Deciding on simulation settings for non-linear models


#### Baseline plot
qqq = 40 
p1 <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) + 
  xlim(0, 10) + theme_bw() + 
  xlab(expression("X"[i]))

##### TOP PANEL
nonlinear_function1 <- function(x, q, gamma1) {
  1-gamma1*x*(x-10)
}
quad_y <- p1 + 
  stat_function(fun = nonlinear_function1, 
                args = list(q = qqq, gamma1 = 0.25), aes(color = "0.25")) + 
  stat_function(fun = nonlinear_function1, 
                args = list(q = qqq, gamma1 = 0.5), aes(color = "0.5")) + 
  stat_function(fun = nonlinear_function1, 
                args = list(q = qqq, gamma1 = 0.125), aes(color = "0.125")) + 
  ylab(expression(mu[i][1])) +
  scale_color_manual(name = expression(paste("Quadratic\ntrend", gamma)), 
                     values =c("0.25" ="black", 
                               "0.5" = "blue", 
                               "0.125" = "gray"))
quad_y
non_linear_function_to_relative_abundance <- function(fun, x, q, ...) {
  fun(x, ...) / (q - 1 + fun(x, ...))
}
quad_relative <- p1 + 
  stat_function(fun = non_linear_function_to_relative_abundance, 
                args = list(fun = nonlinear_function1, 
                            q = qqq,
                            gamma1 = 0.25), aes(color = "0.25"))  + 
  stat_function(fun = non_linear_function_to_relative_abundance, 
                args = list(fun = nonlinear_function1, 
                            q = qqq,
                            gamma1 = 0.5), aes(color = "0.5"))  + 
  stat_function(fun = non_linear_function_to_relative_abundance, 
                args = list(fun = nonlinear_function1, 
                            q = qqq,
                            gamma1 = 0.125), aes(color = "0.125")) + 
  ylab(expression("Z"[i][1])) +
  scale_color_manual(name = expression(paste("Quadratic\ntrend", gamma)), 
                     values =c("0.25" ="black", 
                               "0.5" = "blue", 
                               "0.125" = "gray"))
quad_relative

#### LOWER PANEL
nonlinear_function2 <- function(x, q, gamma1) {
  3 + 0.25*exp(gamma1*x)
}
exp_y <- p1 + 
  stat_function(fun = nonlinear_function2, 
                args = list(q = qqq, gamma1 = 0.25), col = "black") + 
  stat_function(fun = nonlinear_function2, 
                args = list(q = qqq, gamma1 = 0.32), col = "blue") + 
  stat_function(fun = nonlinear_function2, 
                args = list(q = qqq, gamma1 = 0.125), col = "gray") + 
  ylab(expression(mu[i][1])) +
  scale_color_manual(name = expression(paste("Exponential\ntrend", gamma)), 
                     values =c("0.25" ="black", 
                               "0.32" = "blue", 
                               "0.125" = "gray")) 
exp_y
exp_relative <- p1 + 
  stat_function(fun = non_linear_function_to_relative_abundance, 
                args = list(fun = nonlinear_function2, 
                            q = qqq,
                            gamma1 = 0.25), col = "black")  + 
  stat_function(fun = non_linear_function_to_relative_abundance, 
                args = list(fun = nonlinear_function2, 
                            q = qqq,
                            gamma1 = 0.32), col = "blue")  + 
  stat_function(fun = non_linear_function_to_relative_abundance, 
                args = list(fun = nonlinear_function2, 
                            q = qqq,
                            gamma1 = 0.125), col = "gray") + 
  ylab(expression("Z"[i][1])) +
  scale_color_manual(name = expression(paste("Exponential\ntrend", gamma)), 
                     values =c("0.25" ="black", 
                               "0.32" = "blue", 
                               "0.125" = "gray"))
exp_relative

### Shannon 
non_linear_function_to_shannon <- function(fun, x, q, ...) {
  first_relative_abundance <- non_linear_function_to_relative_abundance(fun = fun, x = x, q = q, ...)
  other_relative_abundances <- (1 - first_relative_abundance) / (q - 1)
  
  props <- cbind(first_relative_abundance, 
                 matrix(other_relative_abundances, 
                        ncol = q-1, nrow = length(x), byrow=F))
  apply(X=props, 1, breakaway::true_shannon)
}
exp_shannon <- p1 + 
  stat_function(fun = non_linear_function_to_shannon, 
                args = list(fun = nonlinear_function2, 
                            q = qqq,
                            gamma1 = 0.25), aes(color = "0.25"))  + 
  stat_function(fun = non_linear_function_to_shannon, 
                args = list(fun = nonlinear_function2, 
                            q = qqq,
                            gamma1 = 0.32), aes(color = "0.32"))  + 
  stat_function(fun = non_linear_function_to_shannon, 
                args = list(fun = nonlinear_function2, 
                            q = qqq,
                            gamma1 = 0.125), aes(color = "0.125")) + 
  ylab(expression("Shannon")) +
  scale_color_manual(name = expression(paste("Exponential\ntrend", gamma)), 
                     values =c("0.25" ="black", 
                               "0.32" = "blue", 
                               "0.125" = "gray"))
exp_shannon
quad_shannon <-  p1 +
  stat_function(fun = non_linear_function_to_shannon, 
                args = list(fun = nonlinear_function1, 
                            q = qqq,
                            gamma1 = 0.25), aes(color = "0.25"))  + 
  stat_function(fun = non_linear_function_to_shannon, 
                args = list(fun = nonlinear_function1, 
                            q = qqq,
                            gamma1 = 0.5), aes(color = "0.5"))  + 
  stat_function(fun = non_linear_function_to_shannon, 
                args = list(fun = nonlinear_function1, 
                            q = qqq,
                            gamma1 = 0.125), aes(color = "0.125")) +
  ylab(expression("Shannon")) +
  scale_color_manual(name = expression(paste("Quadratic\ntrend", gamma, sep = "")), 
                     values =c("0.25" ="black", 
                               "0.5" = "blue", 
                               "0.125" = "gray"))
quad_shannon

quad1 <- ggpubr::ggarrange(quad_y, quad_relative, quad_shannon, 
                           nrow = 1, ncol = 3, common.legend=TRUE, legend = "right")
quad2 <- ggpubr::ggarrange(exp_y, exp_relative, exp_shannon, 
                           nrow = 1, ncol = 3, common.legend=TRUE, legend = "right")
shapes <- ggpubr::ggarrange(quad1, 
                            quad2, 
                            nrow = 2)
shapes
ggsave("nonlinear-fig-shapes.pdf", width = 6, height = 3, units="in")
