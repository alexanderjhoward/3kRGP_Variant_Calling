#!/bin/bash -login
#SBATCH --account=adamgrp             # account
#SBATCH -p low2                       # partition
#SBATCH -c 1                          # 1 cores
#SBATCH -t 0:30:00                    # 30 minute reservation
#SBATCH --mem=32000                   # 32 Gb memory

# Initialize conda
. ~/mambaforge/etc/profile.d/conda.sh

# Activate conda environment
conda activate 3kRGP

# Fail on weird errors
set -e
set -x

# Download and index reference genome
wget -P Source/ https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_genome.fasta.gz
gunzip Source/IRGSP-1.0_genome.fasta.gz
samtools faidx Source/IRGSP-1.0_genome.fasta.gz

# Download reference genome annotations
wget -P Source/ https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_2024-07-12.tar.gz
tar -xvzf Source/IRGSP-1.0_representative_2024-07-12.tar.gz
rm Source/IRGSP-1.0_representative_2024-07-12.tar.gz

# Print out final statistics about resource use before job exits
scontrol show job ${SLURM_JOB_ID}
sstat --format 'JobID,MaxRSS,AveCPU' -j ${SLURM_JOB_ID}
