#!/usr/bin/env Rscript
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

df <- read_delim(args[1], col_names=c('a','b'), delim='\t') %>%
    mutate(c = b * 2 )

ggplot(df, aes(a, c)) + geom_point()
ggsave(args[2])