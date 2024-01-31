library(tidyverse)

setwd("C:/Users/David/Downloads")

df <- read.csv("combined_output.csv")

t.test(df %>% filter(autodiff_mode))

lm(time ~ autodiff_mode + exp_type + n + autodiff_mode:exp_type, data=df) %>%
  summary
lm(error ~ autodiff_mode + exp_type + n + autodiff_mode:exp_type, data=df) %>%
  summary

df %>%filter(autodiff_mode == "for") %>% pull(time) %>% mean
df %>%filter(autodiff_mode == "rev") %>% pull(time) %>% mean

df %>%filter(exp_type == "yield") %>% pull(error) %>% median
df %>%filter(exp_type == "trap") %>% pull(error) %>% median
