# Step 2: Call sample genome variants relative to reference genome
The 3kRGP is hosted online and can be queried to pull small genomic regions instead of having to download the entire genome assembly. After pulling these regions, we can compare them to the Nipponbare reference genome and call sequence variants. These variants can then be saved to FASTA files and analyzed down the road.

## Call sequence variants
First, download the list of all 3,024 genome names found in the 3kRGP dataset.

```{bash}
wget -P Source/ https://3kricegenome.s3.amazonaws.com/MANIFEST
grep 'Nipponbare' Source/MANIFEST | awk -F '/' '{print $5}' | sed 's/.realigned.bam//' > Source/samples.txt
```

The **call_3krgp.sh** script does the following:
1. Downloads a BAM genome alignment file from the 3kRGP dataset (just for the genomic regions we defined in Step 1).
2. Calls genomic variants relative to the Nipponbare genome, and filters for high-quality variant calls.
3. Applies the called variants to the reference genome sequence.
4. Saves the mRNA, CDS, and translated peptide sequences for each region of interest.

I prefer to use the UC Davis HPC (FARM) to do this because my account has more storage/processing power than my laptdop does. The **call_3krgp.sh** script is formatted to fit their SLURM queuing system and I input it as an "sbatch" request. I use a loop to automatically input the requests for me.

```bash
while read i; do
	sbatch Scripts/call_3krgp.sh $i
done < Source/samples.txt
```

The variant call files (.vcf) are saved under the "Output/VCF" directory. The genomic files (.fa) are saved under "Output/FASTA" directory. The peptide files (.pep) are saved under "Output/PEP".
