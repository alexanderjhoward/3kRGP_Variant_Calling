# Step 1: Identify regions of interest in reference genome

The 3kRGP has sample genomes aligned against the Nipponbare IRGSP-1.0 reference genome. We need to query these alignments to get sequence variants. Since the sample genomes were assembled relative to Nipponbare, we first need to figure out where in the Nipponbare genome our genes of interest are located. These will be the locations that we then retrieve from each sample genome down the line.

## Download the Nipponbare genome assembly and annotation
First, we need the reference genome and annotations.

```{bash}

	sbatch Scripts/download_reference.sh

```

## Set up the IGV genome browser
Once everything is downloaded, set up an IGV profile to look at sample alignments. [IGV](https://igv.org/doc/desktop/#DownloadPage/) is a genome browser that helps visualize genome annotations and sequence alignments.
1. Open IGV and select "Genomes > Load genome from file”. Select the **IRGSP-1.0_genome.fa** file (downloaded to the "Source" directory).  Once loaded on IGV, the IRGSP-1.0 genome should show up in the browser and be available to select from the drop-down list of genomes.
2. Load in genome annotations by selecting “File > Load from file” and choosing your annotation file. I like the **transcripts.gff** file (found in the "Source/IRGSP-1.0_representative" directory).
3. After setting your browser up, you can save your session by creating an IGV profile. Select “File > Save session” and name the session file whatever you’d like. Save this file to a directory you will remember (such as the "Source" directory). If you add any new annotations or alignments to this session, make sure to re-save the session file before closing IGV! To re-load your session, go to “File > Open session” and select the .xml session file you saved.

*Note: If you move your genome, annotations, sequence files, or .xml session file to a different directory then the session will no longer load properly in IGV, so make sure these files stay where they are when you save your session.*

## Locate genes of interest in the Nipponbare genome
Next, use BLAST to search for your genes of interest within the Nipponbare reference genome. Input all your sequences into with script with a single FASTA file (Mine is called **OsPSY_vars.fa**, found in the "Source" directory. You can use any other FASTA you prefer).
```{bash}

	sbatch Scripts/blast_reference.sh Source/OsPSY_vars.fa

```

The **blast_cleanup.R** script is used next to find the top BLAST hit for each sequence of interest. These are output as a spreadsheet called **IRGSP-1.0_IGVlocs.csv** in the "Output" directory.

```{r}

	Scripts/blast_cleanup.R

```

<center>
<img src="Output/Figures/Before_Annotation.png">
</center>

With this set of genomic regions, we should check that they look correct. My genes of interest came from Kitaake, so I aligned my top BLAST hits against the Kitaake sequences to make sure they looked very similar. I also used IGV to see what the transcript name for each gene region was. Below are my results.

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

I added these transcript names to the **IRGSP-1.0_IGVlocs.csv** spreadsheet by their corresponding sequence.

<center>
<img src="Output/Figures/After_Annotation.png">
</center>
