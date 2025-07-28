This directory contains all related to structural variant (SV) calling. 
Note that because of the significanly better chromosomal resolution, for SV calling in the HiFi data, the T2T reference only was used throughout.
The tools used were NanomonSV, SAVANA, Sawfish and Severus, which are each organised by a separate subdirectory in this folder which include:
the .yaml files for each conda environment (each caller has a separate conda environment) and the scripts for each. Also note that the 
scripts for somatic callers have the same setup metrics. 

Since Sawfish and Severus can be run also on the unpaired sample, there are scripts for structural variant discovery on the hg38 reference. Note that this 
was only done to demonstrate that the CIC breakpoints of interest are NOT discovered when calling on the hg38 reference (e.g. the ChrUn is called instead of chr1). 

There are a couple of caveats connected to pipelining here: 
1) nanomonsv output is not sorted so before annotation it has to be sorted with bcfgools using the post-caller_filtering.sh script
2) sawfish is not a somatic caller so intersection by bcftools has to be performed, as well as normalisation, sorting and indexing
3) sometimes this fails so some incorrect variants where the start position is bigger than the end (e.g. inversions) have to be excluded - does not happen always
4) severus outputs a tsv file (ugly and inappropriate), which has to first be converted to a vcf before annotation - there is a helper function for this in the variant_callers_filtering.ipynb script
5) sawfish ouput is also used for germline variant annotation - e.g. the blood sample is good for this since there shouldn't be somatic variants which are in the tumor
6) annotation is carried out by firsly converting all the Liftoff Refseq annotations gff3 file (manually reviewed genic locations) to a bed file which is used for annotation
7) bcftools annotation can only do one file at a time so first, annotation with the RefSeq and then with RepeatMasker is carried out for each of the vcf file
8) For the ease of running downstream filtering, I have moved all the bcftools annotation output files to one folder called base_dir = "/lustre/scratch125/cellgen/behjati/lr26/SV/" in the downstream Python variant_callers_filtering.ipynb script
9) Then, all the other filtering functions apart from the helper for Severus conversion to vcf can be run as downstream analysis - in the variant_callers_filtering.ipynb script

All the outputs of the variant_callers_filtering.ipynb are included as .csv files or alternatively .xlsx files in the respective directories. Consensus calls are for INS, DEL, BND are also included.
