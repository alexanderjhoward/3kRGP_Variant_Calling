#!/bin/bash -login
#SBATCH --account=adamgrp             # account
#SBATCH -p high2                      # partition
#SBATCH -c 1                          # 1 cores
#SBATCH -t 3:00:00                    # 3 hour reservation
#SBATCH --mem=32000                   # 32 Gb memory

# Initialize conda
. ~/mambaforge/etc/profile.d/conda.sh

# Actiate conda environment
conda activate 3kRGP

# Fail on weird errors
set -e
set -x

# Define input variables
BAM=${1}

# Download sample sequence realignment against Nipponbare genome
wget https://3kricegenome.s3.us-east-1.amazonaws.com/Nipponbare/${BAM}.realigned.bam
wget https://3kricegenome.s3.us-east-1.amazonaws.com/Nipponbare/${BAM}.realigned.bam.bai

while read a b c d; do

	# Subset out genomic regions of interest from alignment file
	samtools view -b -h ${BAM}.realigned.bam ${a}:${b}-${c} > ${BAM}_${d}.bam

done < ../1_Reference_Genome/Output/OsPSY_locs.txt

rm ${BAM}.realigned.bam
rm ${BAM}.realigned.bam.bai

while read a b c d; do
	
	# Call variants with FreeBayes
	freebayes -f ../1_Reference_Genome/Source/IRGSP-1.0_genome.fasta ${BAM}_${d}.bam > Output/VCF/${BAM}_${d}.vcf
	rm ${BAM}_${d}.bam

	# Filter low quality calls
        vcftools --vcf Output/VCF/${BAM}_${d}.vcf --minQ 20 --recode --recode-INFO-all --out Output/VCF/${BAM}_${d}_q20
        rm Output/VCF/${BAM}_${d}.vcf

	# Compress and index VCF
	bgzip Output/VCF/${BAM}_${d}_q20.recode.vcf
	tabix -p vcf Output/VCF/${BAM}_${d}_q20.recode.vcf.gz

	# Apply variants to region of interest
	samtools faidx ../1_Reference_Genome/Source/IRGSP-1.0_genome.fasta ${a}:${b}-${c} | bcftools consensus Output/VCF/${BAM}_${d}_q20.recode.vcf.gz -o Output/FASTA/${BAM}_${d}.fa

	# Edit FASTA header to give each sequence a unique name
	sed -i "1 s/^.*$/>${BAM}_${d}/" Output/FASTA/${BAM}_${d}.fa

done < ../1_Reference_Genome/Output/OsPSY_locs.txt

# Print out final statistics about resource use before job exits
scontrol show job ${SLURM_JOB_ID}
sstat --format 'JobID,MaxRSS,AveCPU' -j ${SLURM_JOB_ID}
