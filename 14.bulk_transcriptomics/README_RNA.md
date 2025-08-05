This directory contains all bulk transcriptomic scripts, including the feature table used in tximport (R). 
There are separate STAR and Salmon scripts for our samples, which involved re-extracting reads from a cram file. 
Note that the extraction is in the star-our-patient.sh, which means that the salmon-our-samples.sh should be run after 
to use the extracted fastq files. Similarly, fastqc-fastq.gz.sh is a script to run fastqc on all the samples obtained. 
Note that extrernal samples were obtained in either the zipped or unzipped fastq form, so that salmon-cohort.sh and star-cohort-fin.sh 
can be run rightaway without the need to re-extract. 

Also note that bulk-rna-blats.sh include the old pathways, as this process was performed with the old pathway system and psl files were viewed for matches. 

There is a downstream RNA-seq-deseq2.ipynb script which is written in R and should be run after all Salmon quant ouputs were generated. 
There is a csv file which includes significant genes and pathways found in the downstream processing of gene expression. 
