library(tidyverse)

#aeds <- read_tsv('output/summary/aeds.tsv')
aeds <- read_tsv(snakemake@input[[1]])

aeds %>% 
  ggplot(aes(x=AED, y=cumulative_freq))+
  geom_line(aes(colour=Strain)) +
  xlab('Annotation Edit Distance') +
  ylab('Cumulative number of transcripts')

ggsave(snakemake@output[[1]])