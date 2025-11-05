#Setup
library(rstudioapi)
library(tidyverse)
setwd(dirname(getActiveDocumentContext()$path))

# Load in BLAST results
df <- read.csv('../Output/OsPSY_search.txt', header = F)
colnames(df) <- c('qseqid','sseqid','evalue','bitscore','qstart','qend','qseq','sstart','send','sseq')

# Filter for sequences of interest
filtered <- df %>%
  group_by(qseqid) %>%
  filter(evalue <= min(evalue) & qstart == 1) %>%
  mutate(strand = ifelse(sstart < send, "+", "-")) %>%
  mutate(start = ifelse(sstart < send, sstart, send)) %>%
  mutate(end = ifelse(sstart < send, send, sstart)) %>%
  mutate(igv = paste0(sseqid, ":", start, "-", end))

# Format output file to search for corresponding transcript ID in IGV
output <- filtered %>% select(qseqid, igv)
write.csv(output, file="../Output/IRGSP-1.0_OsPSY_IGVlocs.csv", quote=FALSE, row.names = FALSE)
