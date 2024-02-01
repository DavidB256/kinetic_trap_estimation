# Setup
library(tidyverse)
setwd("C:/Users/David/Desktop/Johnson/R")

# Import data from CSV output by `reformat_output.ipynb`
df <- read.csv("combined_output.csv")

# Exploratory data analysis
df %>%filter(autodiff_mode == "for") %>% pull(time) %>% mean
df %>%filter(autodiff_mode == "rev") %>% pull(time) %>% mean

df %>%filter(exp_type == "yield") %>% pull(error) %>% median
df %>%filter(exp_type == "trap") %>% pull(error) %>% median

# Test for effects of autodiff mode (forward vs. reverse) and experimental data 
# type (yield vs. trap) on runtime while controlling for reaction network size (`n`)
# and the interaction of mode and data type.
lm(time ~ autodiff_mode + exp_type + n + autodiff_mode:exp_type, data=df) %>%
  summary
# Ditto, but for mean-squared error in estimates instead of runtime.
lm(error ~ autodiff_mode + exp_type + n + autodiff_mode:exp_type, data=df) %>%
  summary
