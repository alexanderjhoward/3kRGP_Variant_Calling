# Step 1: Identify regions of interest in reference genome

The 3kRGP provides sample genomes aligned against the Nipponbare IRGSP-1.0 reference genome. These alignments are what we query to get sequence variants. Since the sample genomes were assembled relative to Nipponbare, we need to figure out where in the Nipponbare genome our genes of interest are located. These are the locations that will then be retrieved from each sample genome down the line.

First, we download the reference genome and annotations.

```{bash}

	sbatch Scripts/download_reference.sh

```

Once downloaded, we can set up an IGV profile to look at sample alignments later on. [IGV](https://igv.org/doc/desktop/#DownloadPage/) is a genome browser that helps visualize genome annotations and sequence alignments.
1. Open IGV and select "Genomes” > “Load genome from file”. Select the **IRGSP-1.0_genome.fa** file (which should have been downloaded to the "Source" directory).  Once loaded into IGV, the IRGSP-1.0 genome should appear in the browser and be available to select in the browser drop-down list of genomes.
2. Load in any desired annotations by selecting “File” > “Load from file” and choosing your annotation file. I like the **transcripts.gff** file (found in the "Source/IRGSP-1.0_representative" directory).
3. After setting this up, you can save your work by creating an IGV profile. Select “File” > “Save session” and name the session file whatever you’d like and save it to a directory you will remember (such as the "Source" directory). If you add any new annotations or alignments to this session, be sure to re-save the session file before closing IGV! To re-load your session, go to “File” > “Open session” and select the .xml session file you saved. Note that if you move your genome/annotations/sequence files or the .xml file to different directories after saving then it will no longer load properly, so make sure these files all stay put where they are from now on.

Next, we use BLAST to search for our genes of interest within the Nipponbare reference genome. Make sure to input all your sequences as a single FASTA file (I have mine in **OsPSY_vars.fa** under the "Source" directory but you can provide any file you prefer).
```{bash}

	sbatch Scripts/blast_reference.sh Source/OsPSY_vars.fa

```

Finally, we can investigate the top alignment matches returned by BLAST and save the top genomic regions for use in downstream scripts. I ran the **blast_cleanup.R** script to do this.

```{r}

	Scripts/blast_cleanup.R

```

It's a good idea to check the regions we saved and how they compare to our original genes of interest. My genes of interest came from Kitaake, so I aligned the top BLAST hits from Nipponbare against the original Kitaake sequences to make sure they look very similar. I also used IGV to see what the transcript name for each gene region was. Below are my results.

|Sequence|Location|Transcript|% Identity to Kitaake|Strand|
|:---:|:---:|:---:|:---:|:---:|
|OsPSY1|chr05:23958390-23959669|Os05t0487100-01|99.8%|+|
|OsPSY2|chr05:23978454-23979760|Os05t0487300-01|99.8%|+|
|OsPSY3|chr01:9711570-9712578|Os01t0276900-01|98.2%|-|
|OsPSY4|chr01:34674485-34676267|Os01t0815400-01|99.8%|-|
|OsPSY5|chr11:23063010-23064583|Os11t0600600-01|99.9%|+|
|OsPSY6|chr07:26188163-26189948|Os07t0631300-01|99.9%|-|
|OsPSY7|chr05:26922134-26923393|Os05t0542300-01|99.8%|-|
|OsPSY8|chr01:8986139-8987609|Os01t0264400-01|99.4%|+|
