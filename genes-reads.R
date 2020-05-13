library(tidyverse)
summary <- read_csv('output/summary/total.csv')
summary %>% 
  ggplot(mapping = aes(x=Reads, y=Genes)) +
  geom_line(aes(colour=Strain)) + 
  ylim(0, 20000) +
  xlim(0, 10E7) +
  xlab('RNA-seq reads used') +
  ylab('Genes annotated by MAKER')
  #geom_text(hjust=0, vjust=0, nudge_x = 0.02E7)
ggsave(snakemake@output[[1]])

cor.test(summary$Reads, summary$Genes)
