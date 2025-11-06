#!/bin/bash
# Define variables
FASTA=${1}

# Check if variable is empty or unset
if [[ -z "$FASTA" ]]; then
  echo "Error: FASTA is required but not provided." >&2
  exit 1
fi

# Proceed if variable is provided
echo BLAST is querying with "$FASTA"

# BLAST search for genes of interest against Nipponbare reference genome
cp Source/IRGSP-1.0_genome.fasta Source/BLAST/
makeblastdb -in Source/BLAST/IRGSP-1.0_genome.fasta -dbtype nucl
blastn -query $FASTA -db Source/BLAST/IRGSP-1.0_genome.fasta -outfmt '10 qseqid sseqid evalue bitscore qstart qend qseq sstart send sseq' -out Output/IRGSP-1.0_BLASTsearch.txt
