#!/usr/bin/env Rscript

library(tidyverse)

mu <- 1.4e-9
gen <- 5
df.msmc <- read.table(snakemake@input[[1]], header=TRUE)

full.filename <- name <- basename(snakemake@input[[1]])

plot.title <- gsub("_msmc2.final.txt","",full.filename)

df.msmc %>% 
  mutate(YearsAgo = (left_time_boundary/mu)*gen,
         Ne = (1/lambda)/(2*mu)) %>% 
  ggplot(aes(x = YearsAgo, y = Ne)) +
  geom_step() +
  scale_x_log10() +
  theme_bw() + ggtitle(plot.title)

ggsave(snakemake@output[[1]], width = 8.5, height = 8.5)