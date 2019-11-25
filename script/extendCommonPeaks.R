#!/usr/bin/env Rscript
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

h <- paste0(c("chrom", "start", "end", "name", "signal", "strand", 
              "fc", "logp", "logq", "summit"), c(rep("1", 10), rep("2", 10)))
x <- read_delim(args[1], col_names=c(h, 'overlap'), delim='\t')

# extend the overlap
x <- mutate(x, chrom = chrom1,
            start = if_else(start1 < start2, start1, start2),
            end = if_else(end2 > end1, end2, end1),
            name = paste(name1, name2, sep='-'),
            signal = ceiling((signal1 + signal2)/2), strand = '.',
            fc = (fc1 + fc2)/2, 
            logp = (logp1 + logp2)/2,
            logq = (logq1 + logq2)/2, 
            summit = (summit1 + summit2)/2 )

select(x, chrom, start, end, name, signal, strand, fc, logp, logq, summit) %>%
  write_delim(args[2], delim="\t", col_names = FALSE)
