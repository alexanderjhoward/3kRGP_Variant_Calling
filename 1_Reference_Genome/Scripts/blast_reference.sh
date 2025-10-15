#!/bin/bash -login
#SBATCH --account=adamgrp             # account
#SBATCH -p low2                       # partition
#SBATCH -c 1                          # 1 cores
#SBATCH -t 0:15:00                    # 15 minute reservation
#SBATCH --mem=8000                    # 8 Gb memory

# Initialize conda
. ~/mambaforge/etc/profile.d/conda.sh

# Activate conda environment
conda activate 3kRGP

# Fail on weird errors
set -e
set -x

# BLAST search for genes of interest against Nipponbare reference genome
cp Source/IRGSP-1.0_genome.fasta Source/BLAST/
makeblastdb -in Source/BLAST/IRGSP-1.0_genome.fasta -dbtype nucl
blastn -query Source/OsPSY_vars.fa -db Source/BLAST/IRGSP-1.0_genome.fasta -outfmt '10 qseqid sseqid evalue bitscore qstart qend qseq sstart send sseq' -out Output/OsPSY_search.txt

# Print out final statistics about resource use before job exits
scontrol show job ${SLURM_JOB_ID}
sstat --format 'JobID,MaxRSS,AveCPU' -j ${SLURM_JOB_ID}
