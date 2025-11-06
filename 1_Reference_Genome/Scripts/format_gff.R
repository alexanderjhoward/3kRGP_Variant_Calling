#Setup
library(tidyverse)
library(ape)

# Read in genome annotation
gff <- read.gff('Output/genes_of_interest.gff', na.strings = c('.', '?'), GFF3 = TRUE)

# Load in IGV data (after we edited it to include our transcripts of interest)
locs <- read.table('Output/IRGSP-1.0_IGVlocs.txt', col.names=c('qseqid', 'igv', 'transcript'))
transcripts <- locs %>% pull(transcript)

# Pull genes of interest from genome annotation
genes <- gff %>%
  filter(grepl(paste(transcripts, collapse='|'), attributes)) %>%
  mutate(attributes = gsub(";.*", "", attributes)) %>%
  mutate(attributes = gsub(".*=", "", attributes)) %>%
  merge(locs, by.x="attributes", by.y="transcript") %>%
  filter(type %in% c('mRNA')) %>%
  mutate(strand = as.character(strand)) %>%
  mutate(score = "0") %>%
  arrange(qseqid) %>%
  select(seqid, start, end, qseqid, score, strand)
write.table(genes, file='Output/genes.bed', quote=FALSE, sep="\t", row.names = F, col.names = F)

# Pull CDS of interest from genome annotation and order the exons by strandedness
#   - If the strand is forward (+), order exons in ascending order
#   - If the strand if reverse (-), order exons in descending order
pull_cds <- function(a, b){
  locs_subset <- gff %>%
    filter(grepl(a, attributes)) %>%
    mutate(attributes = gsub(";.*", "", attributes)) %>%
    mutate(attributes = gsub(".*=", "", attributes)) %>%
    merge(locs, by.x="attributes", by.y="transcript") %>%
    filter(type %in% c('CDS')) %>%
    mutate(strand = as.character(strand)) %>%
    select(qseqid, type, seqid, start, end, strand)
  strand <- unique(locs_subset %>% pull(strand))
  if(strand=="+"){
    locs_subset <- locs_subset %>% arrange(start)
  } else {
    locs_subset <- locs_subset %>% arrange(desc(start))
  }
  locs_subset <- locs_subset %>%
    mutate(coords = paste0(seqid, ":", start, "-", end)) %>%
    select(coords)
  write.table(locs_subset, file=paste0('Output/', b,'_CDS.txt'), quote=FALSE, sep="\t", row.names = F, col.names = F)
}

for (i in 1:nrow(locs)){
  pull_cds(locs[i,3], locs[i,1])
}
