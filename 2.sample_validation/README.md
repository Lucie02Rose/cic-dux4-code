This section contains scripts related to sample validation. 
Note that firstly, the Allele Integrator (genotypeCheck.R utilised by the wrapper) were run provided that all 
dependencies and paths are stisfied and no files are missing. 
Secondly, a bed file of the sampled locations is created n the allele_integrator_hifi_notebook. 
This bed file is then utilised by the samtools-sn-pb.sh script which samples these locations with bcftools mpileup.

The downstream processing is carried out in the allele_integrator_hifi_notebook, namely file conversion and pariwise heatmap construction.
