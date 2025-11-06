#!/bin/bash
# Download and index reference genome
wget -P Source/ https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_genome.fasta.gz
gunzip Source/IRGSP-1.0_genome.fasta.gz
samtools faidx Source/IRGSP-1.0_genome.fasta

# Download reference genome annotations
wget -P Source/ https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_2024-07-12.tar.gz
tar -xvzf Source/IRGSP-1.0_representative_2024-07-12.tar.gz -C Source/
rm Source/IRGSP-1.0_representative_2024-07-12.tar.gz
