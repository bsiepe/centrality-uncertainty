library(tidyverse)
data_frame(Cauchy = rcauchy(1e6, 0, .1),
           Students_t_2 = rt(1e6, 2) * .1,
           Students_t_3 = rt(1e6, 3) * .1,
           Normal = rnorm(1e6) * .1) %>%
  pivot_longer(cols = everything()) %>%
  ggplot(aes(value, col = name)) +
  geom_density(bins = 1e3, ) +
  xlim(-1, 1) +
  ylim(0,2) +
  theme_minimal() +
  labs(title = "Histogram of 10,000 Cauchy(0,0.1) random variables", x =
         "Value", y = "Frequency")
