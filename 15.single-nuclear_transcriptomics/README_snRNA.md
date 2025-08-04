This directory contains all single-nuclear transcriptomics scripts and results. 
There is a sc_fusions.ipynb script which is used to extract soft clipped sequences from bam files for all 
channels and then these fasta are used in the blat-sc.sh script to find whether breakpoints at the 21th exon 
of the CIC gene correspond. During the BLAT, cell names of the sequences are retained to then use to plot
the fusion positive cells in a UMAP in the snRNA-final.ipynb. This script contains the entire snRNA workflow. 

There are also a number of output files in this directory, including the Leiden clusters, differentially expressed 
genes in the tumor, GSEApy results and all the extracted soft clips in the fasta files (per channel).
