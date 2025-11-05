#Setup
library(rstudioapi)
library(tidyverse)
library(ape)
setwd(dirname(getActiveDocumentContext()$path))

# Read in genome annotation
gff <- read.gff('../Source/IRGSP-1.0_representative/transcripts.gff', na.strings = c(".", "?"), GFF3 = TRUE)

# Load in IGV data (after we edited it to include our transcripts of interest)
locs <- read.csv('../Output/IRGSP-1.0_OsPSY_IGVlocs.csv')
transcripts <- locs %>% pull(transcript)

# NOTE: May want to make this output as a BED format for ease of use downstream.


# Pull transcripts of interest from genome annotation
genes <- gff %>%
  filter(grepl(paste(transcripts, collapse='|'), attributes)) %>%
  mutate(attributes = gsub(";.*", "", attributes)) %>%
  mutate(attributes = gsub(".*=", "", attributes)) %>%
  merge(locs, by.x="attributes", by.y="transcript") %>%
  filter(type %in% c('mRNA', 'CDS')) %>%
  mutate(strand = as.character(strand)) %>%
  arrange(qseqid) %>%
  select(seqid, start, end, qseqid, type, strand)
write.table(genes, file='../Output/test_locs.tsv', quote=FALSE, sep="\t", row.names = F, col.names = F)

# Pull CDSs of interest from genome annotation
# We want to order the exons how they will be translated
# This means if the strand is forward (+) we order exons from lowest bp start to highest
# If the strand if reverse (-) we order exonds fom highest bp start to lowest
pull_cds <- function(a){
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
  return(locs_subset)
}
cds <- NULL
for (i in transcripts){
  out <- pull_cds(i)
  cds <- rbind(cds, out)
}
rm(out, i)

# Note! Sometimes genome annotations can be wrong.
# The IRGSP-1.0 annotation for OsPSY5 is incorrect when consulting the Kitaake annotation.
# Remove the wrong OsPSY5 CDS
cds <- cds %>%
  filter(!(qseqid %in% c('OsPSY5')))

# We can manually define the annotation.
# OsPSY5 Exon1: chr11:23063675-23063771
# OsPSY5 Exon2: chr11:23063896-23063954
# OsPSY5 Exon3: chr11:23064191-23064322
# Note that OsPSY5 is forward (+) so we order exons from lowest to highest bp starts
ospsy5_e1 <- c('OsPSY5', 'Os11t0600600-01', 'CDS', 'chr11', '23063675', '23063771', '+')
ospsy5_e2 <- c('OsPSY5', 'Os11t0600600-01', 'CDS', 'chr11', '23063896', '23063954', '+')
ospsy5_e3 <- c('OsPSY5', 'Os11t0600600-01', 'CDS', 'chr11', '23064191', '23064322', '+')
cds <- rbind(cds, ospsy5_e1, ospsy5_e2, ospsy5_e3)
rm(ospsy5_e1, ospsy5_e2, ospsy5_e3)
