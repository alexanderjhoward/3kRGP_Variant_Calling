# Step 2: Call sample genome variants relative to reference genome
The 3kRGP is hosted online and can be queried to pull small genomic regions instead of having to download the entire genome assembly. After downloading these alignments, we can compare them to the Nipponbare reference genome and call sequence variants. These variants can then be saved to FASTA files and analyzed down the road.

## Call sequence variants
The **call_variants.sh** script does the following:
1. Download a BAM genome alignment file from the 3kRGP dataset according to the list of samples we made (samples.txt)
2. Pull the genomic regions that correspond to our sequences of interest (found in **OsPSY_locs.txt**, generated during Step 1)
3. Call variants of the genome relative to the Nipponbare genome, and filter for high-quality variant calls
4. Apply the variants called to the reference genome sequence and save the sequence as a FASTA file



First we need to download the list of all 3,024 genome names found in the 3kRGP dataset.

```{bash}

	wget -P Source/ https://3kricegenome.s3.amazonaws.com/MANIFEST
	grep 'Nipponbare' Source/MANIFEST | awk -F '/' '{print $5}' | sed 's/.realigned.bam//' > Source/samples.txt

```

Now, with all our genomic locations determined, we can go through the process of downloading sample genomes, saving genomic regions of interest, and calling variants of the sample genome relative to the Nipponbare genome. This is all done through the code shown below. 

I use the “screen” function in Unix to input SLURM requests because I can then run this code in the background without having to keep tabs on it. I set a limit on the loop to input one SLURM request per 30 seconds. This was done to reduce file storage demand and improve download speed. That means 120 requests/hr, so for 3,024 samples this requires 25.2 hours of runtime to complete. I needed hundreds of GB of file storage to accomodate this process. If you have a lot of file storage, this delay in requests could probably be shortened. Alternatively, the entire set of resequenced alignments could be downloaded and stored locally to search against as a database. This would require a few TB of storage, but would make each variant search much faster.

```{bash}

	screen -S PSY_run # Set a session name for the screen to help distinguish it later

	while read i; do
		sbatch Scripts/call_variants.sh $i
		sleep 30s
	done < ../1_Reference_Genome/Output/OsPSY_locs.txt

	## press Ctrl + a + d to "detach" from the screen and let it run in the background

	screen -ls # lists all active screens, you should see one named (string of numbers).PSY_run 
	screen -r 2771512 # Input the string of numbers for the screen session to return to this screen. Make sure to detatch from the screen again afterwards.

```

The **call_variants.sh** script does the following:
1. Download a BAM genome alignment file from the 3kRGP dataset according to the list of samples we made (samples.txt)
2. Pull the genomic regions that correspond to our sequences of interest (found in **OsPSY_locs.txt**, generated during Step 1)
3. Call variants of the genome relative to the Nipponbare genome, and filter for high-quality variant calls
4. Apply the variants called to the reference genome sequence and save the sequence as a FASTA file

The variant call files (VCFs) are saved under the "Output/VCF" directory. The genomic files (FASTAs) are saved under "Output/FASTA" directory.
