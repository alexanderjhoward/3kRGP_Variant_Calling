#Setup
library(dplyr)

# Load in BLAST results
df <- read.csv('Output/IRGSP-1.0_BLASTsearch.txt', header = F)
colnames(df) <- c('qseqid','sseqid','evalue','bitscore','qstart','qend','qseq','sstart','send','sseq')

# Filter for sequences of interest
filtered <- df %>%
  group_by(qseqid) %>%
  filter(evalue <= min(evalue) & qstart == 1) %>%
  mutate(strand = ifelse(sstart < send, "+", "-")) %>%
  mutate(start = ifelse(sstart < send, sstart, send)) %>%
  mutate(end = ifelse(sstart < send, send, sstart)) %>%
  mutate(igv = paste0(sseqid, ":", start, "-", end)) %>%
  mutate(transcript=NA)

# Output CSV table with IGV locations of genes of interest and a blank space for what transcript corresponds
output <- filtered %>% select(qseqid, igv, transcript)
write.table(output, file="Output/IRGSP-1.0_IGVlocs.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
