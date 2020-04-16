if(!require(tidyverse)){
  install.packages('tidyverse')
  library(tidyverse)
}
read_gff <- function(filename){
  gff = read_tsv(filename, col_names=F, comment='#') %>% 
    rename(seqid = X1,
           source = X2,
           type = X3,
           start = X4,
           end = X5,
           score = X6,
           strand = X7,
           phase = X8,
           attributes = X9) %>% 
    mutate(seqid = as_factor(seqid),
           source = as_factor(source),
           type = as_factor(type),
           strand = as_factor(strand),
           phase = as_factor(phase))
  return(gff)
}
maker_output <- read_gff(snakemake@input[[1]])
#reference_annotations <- read_gff('output/ref/annotation/ref-maker_shared_classes.gff')
maker_summary <- maker_output %>% 
  filter(type != 'contig') %>% 
  group_by(type, source) %>%
  summarise(Count = n()) %>% 
  arrange(desc(Count))
print(maker_summary)
maker_summary %>%
  rename(Type=type, Evidence=source) %>% 
  ggplot(aes(x=Type, y=Count)) +
  geom_col(aes(fill=Evidence)) +
  ggtitle(paste(snakemake@wildcards[['strain']], ' feature types and evidences'))
ggsave(snakemake@output[[1]], device='png')

#reference_summary <- reference_annotations %>% 
#  group_by(type, source) %>% 
#  summarise(n = n()) %>% 
#  arrange(n)
#reference_summary %>% 
#  ggplot(aes(x=type, y=n)) + geom_col(aes(fill=source))

#maker_output %>% 
#  filter(type != 'contig') %>% 
#  ggplot(aes(x = type)) + geom_bar(aes(fill = source))
#reference_annotations %>% 
#  filter(type != 'contig') %>% 
#  ggplot(aes(x = type)) + geom_bar(aes(fill = source))