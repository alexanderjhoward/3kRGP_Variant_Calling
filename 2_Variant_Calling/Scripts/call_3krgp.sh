#!/bin/bash -login
#SBATCH --account=adamgrp             # account
#SBATCH -p low2                       # partition
#SBATCH -c 1                          # 1 cores
#SBATCH -t 0:10:00                    # 10 minute reservation
#SBATCH --mem=4000                    # 4 Gb memory

# Initialize and activate conda environment
. ~/mambaforge/etc/profile.d/conda.sh
conda activate 3kRGP

# Fail on weird errors
set -e
set -x

# Define input variables
BAM=${1}

# Subset out regions of interest from alignment file (aligned relative to Nipponbare)
samtools view -b -h http://s3.amazonaws.com/3kricegenome/Nipponbare/${BAM}.realigned.bam -M -L ../1_Reference_Genome/Output/genes.bed > Output/BAM/${BAM}_subset.bam
rm ${BAM}.realigned.bam.bai

# Call variants relative to Nipponbare genome with FreeBayes
freebayes -f ../1_Reference_Genome/Source/IRGSP-1.0_genome.fasta Output/BAM/${BAM}_subset.bam > Output/VCF/${BAM}_subset.vcf
rm Output/BAM/${BAM}_subset.bam

# Filter out low quality calls
vcftools --vcf Output/VCF/${BAM}_subset.vcf --minQ 20 --recode --recode-INFO-all --out Output/VCF/${BAM}_subset_q20
rm Output/VCF/${BAM}_subset.vcf

# Compress and index VCF
bgzip Output/VCF/${BAM}_subset_q20.recode.vcf
tabix -p vcf Output/VCF/${BAM}_subset_q20.recode.vcf.gz

# Call variants within regions of interest
while read a b c d e f; do

	# Reverse-stranded variant calling
	if [[ $f == "-" ]]; then
		# mRNA
                samtools faidx ../1_Reference_Genome/Source/IRGSP-1.0_genome.fasta ${a}:${b}-${c} | bcftools consensus Output/VCF/${BAM}_subset_q20.recode.vcf.gz | seqtk seq -r > Output/FASTA/${BAM}_${d}.fa
		sed -i "1 s/^.*$/>${BAM}_${d}_mRNA/" Output/FASTA/${BAM}_${d}.fa
		# CDS
                samtools faidx ../1_Reference_Genome/Source/IRGSP-1.0_genome.fasta -r ../1_Reference_Genome/Output/${d}_CDS.txt | bcftools consensus Output/VCF/${BAM}_subset_q20.recode.vcf.gz | seqtk seq -r | sed '1!{/>/d}' > Output/FASTA/${BAM}_${d}_CDS.fa
                sed -i "1 s/^.*$/>${BAM}_${d}_CDS/" Output/FASTA/${BAM}_${d}_CDS.fa
		# Protein
                transeq Output/FASTA/${BAM}_${d}_CDS.fa Output/PEP/${BAM}_${d}.pep
		sed -i "1 s/^.*$/>${BAM}_${d}_PEP/" Output/PEP/${BAM}_${d}.pep
	# Forward-stranded variant calling
	elif [[ $f == "+" ]]; then
		# mRNA
		samtools faidx ../1_Reference_Genome/Source/IRGSP-1.0_genome.fasta ${a}:${b}-${c} | bcftools consensus Output/VCF/${BAM}_subset_q20.recode.vcf.gz -o Output/FASTA/${BAM}_${d}.fa
		sed -i "1 s/^.*$/>${BAM}_${d}_mRNA/" Output/FASTA/${BAM}_${d}.fa
		# CDS
		samtools faidx ../1_Reference_Genome/Source/IRGSP-1.0_genome.fasta -r ../1_Reference_Genome/Output/${d}_CDS.txt | bcftools consensus Output/VCF/${BAM}_subset_q20.recode.vcf.gz | sed '1!{/>/d}' > Output/FASTA/${BAM}_${d}_CDS.fa
		sed -i "1 s/^.*$/>${BAM}_${d}_CDS/" Output/FASTA/${BAM}_${d}_CDS.fa
		# Protein
		transeq Output/FASTA/${BAM}_${d}_CDS.fa Output/PEP/${BAM}_${d}.pep
		sed -i "1 s/^.*$/>${BAM}_${d}_PEP/" Output/PEP/${BAM}_${d}.pep
	else
		echo "ERROR: No strandedness provided"
		break
	fi

done < ../1_Reference_Genome/Output/genes.bed

# Collect all sample results into single files
cat Output/FASTA/${BAM}_*.fa > Output/FASTA/${BAM}.fa
rm Output/FASTA/${BAM}_*.fa
cat Output/PEP/${BAM}_*.pep > Output/PEP/${BAM}.pep
rm Output/PEP/${BAM}_*.pep

# Print out final statistics about resource use before job exits
scontrol show job ${SLURM_JOB_ID}
sstat --format 'JobID,MaxRSS,AveCPU' -j ${SLURM_JOB_ID}
