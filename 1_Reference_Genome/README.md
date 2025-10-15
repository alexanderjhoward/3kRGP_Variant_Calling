# Step 1: Identify regions of interest in reference genome

The 3kRGP provides sample genomes aligned against the Nipponbare IRGSP-1.0 reference genome. These alignments are what we query to get sequence variants. Since the sample genomes have been assembled relative to Nipponbare, we first need to figure out where in the Nipponbare genome our genes of interest are located so we can pull that genomic region from each sample.

First, we download the reference genome and annotations.

```{bash}

	sbatch Scripts/download_reference.sh

```

Once downloaded, we can set up an IGV profile to look at sample alignments down the line. [IGV](https://igv.org/doc/desktop/#DownloadPage/) is a genome browser that helps visualize genome annotations and sequence alignments.
1. Open IGV and select "Genomes” > “Load genome from file”. Select the IRGSP-1.0 genome file (which should have been downloaded to the "Source" directory).  Once loaded into IGV, the IRGSP-1.0 genome should appear in the browser and be available to select in the browser drop-down list of genomes.
2. Load in any desired annotations by selecting “File” > “Load from file” and choosing your annotation file. I like the transcripts.gff file (found in the "Source/IRGSP-1.0_representative" directory).
3. With this set up, you can save your work by creating an IGV profile. Select “File” > “Save session” and name the session file whatever you’d like. If you add any new annotations or alignments to the session that you want to load back in later on, be sure to re-save the session before closing IGV! To load your saved session, go to “File” > “Open session” and select the .xml file. Note that if you move your genome/annotations/sequence files or the .xml file to different directories after saving the profile then it will no longer reload properly.

Next, we use BLAST to search for our genes of interest within the Nipponbare reference genome. Make sure you have all your sequences of interest in a single FASTA file (I have mine in "OsPSY_vars.fa" under the "Source" directory).
```{bash}

	sbatch Scripts/blast_reference.sh

```

Finally, we can investigate the top alignment matches returned by BLAST and format the regions for downstream scripts. I ran the **blast_cleanup.R** script to do this.

It's a good idea to check which sequences we're looking at, and how they compare to our reference genes. My reference genes came from Kitaake, so I aligned the top BLAST hits from Nipponbare against the Kitaake sequences to make sure they are very similar. Below are the results I got.

|Sequence|Location|Transcript|% Identity to Kitaake|Strand|
|:---:|:---:|:---:|:---:|:---:|
|OsPSY1|chr05:23958390-23959669|Os05t0487100-01|99.8%|+|
|OsPSY2|chr05:23978454-23979760|Os05t0487300-01|99.8%|+|
