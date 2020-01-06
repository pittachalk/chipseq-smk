#!/usr/bin/env Rscript
library(tidyverse)

df <- read_delim(snakemake@input[[1]], col_names=c('a','b'), delim='\t') %>%
    mutate(c = b * 2 )

ggplot(df, aes(a, c)) + geom_point()
ggsave(snakemake@output[[1]])