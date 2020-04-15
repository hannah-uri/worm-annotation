if(!require(tidyverse)){
  install.packages('tidyverse')
  library(tidyverse)
}

trna <- read_tsv(snakemake@input[[1]])
trna
trna %>% ggplot(aes(x = tRNA, y=fraction, col=sample)) +
  geom_point(size =0.8, alpha = 0.8, position = position_dodge(width = 0.3, preserve = 'total')) +
  scale_y_continuous(trans=snakemake@params[['transformation']]) +
  coord_polar() +
  ylab(paste(snakemake@params[['transformation']], '(fraction of total tRNA expresssion)')) + 
  ggtitle(paste('tRNA Expression (', snakemake@wildcards[['type']], ', r', snakemake@wildcards[['replicate']], ')', sep='')
)
ggsave(snakemake@output[[1]], device='png')