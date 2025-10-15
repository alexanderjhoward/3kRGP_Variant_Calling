# Variant calling with resequencing data from the 3kRGP dataset

This repository outlines my methodology for retreiving all unique genomic variants of a gene of interest from the [3,000 Rice Genomes Project (3kRGP)](https://iric.irri.org/projects/3000-rice-genomes-project) dataset. 

The 3kRGP dataset provides resequencing data (sample reads aligned against a reference genome) and VCF files (variant calls relative to the reference genome). While the VCF file is more convenient to use, I personally have had trouble retreiving indel variants with this data. Calling variants directly from the resequencing data is more computationally taxing, but provides the most control over the entire process.

My method for this approach is as follows:
1. Identify the genomic location of your gene of interest within the Nipponbare reference genome.
2. Download a sample's resequencing data and isolate the genomic region of interest.
3. Call sequence variants within this region relative to the Nipponbare reference genome.
4. Convert the variant calls into DNA and protein sequences.
5. Analyze and visualize the sequence variant results.

