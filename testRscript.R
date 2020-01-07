#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(readr)

df <- read_delim(snakemake@input[[1]], col_names=c('a','b'), delim='\t') %>%
    mutate(c = b * 2 )

pdf(file = snakemake@output[[1]])
ggplot(df, aes(a, c)) + geom_point(colour=snakemake@params[["color"]])
dev.off()

sessionInfo()