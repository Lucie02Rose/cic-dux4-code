# cic-dux4-code
This GitHub repository contains all code (Bash/Python/R), files and outputs (tsv/xlsx/csv/txt/pdf/fasta) from the analyses conducted as part of the CIC-DUX4 fusion sarcoma 
project at The Wellcome Sanger Institute in Professor Sam Behjati's Group.
The repository is structured numerically (1-15) as methods appear in the thesis.
There are yaml files which reslect the conda environments and necessary installations available in each particular section. The overall R, Python and 
base conda environment .yaml files are standalone in this repository.

Throughout the repositories, it will be noticable that home directories differ: e.g. nfs or casm/team274sb to cellgen/behjati. This reflects the migration
caused by the Sanger datacentre incident and directories are stored as variables. 

The Wellcome Sanger Institute uses the LSF system, meaning that all the .sh scripts in all directories would be submitted by using bsub < script. 
All notebooks were worked on using the interactive JupyterHub, with setting up the corresponding Conda Python and R environemnts including installing kernels.
Importantly, all Python and R packages must be installed with either conda (preferrable) or pip with Python, or through conda (preferrable) or compiling in R. 
For R, conda worked for me, including bioconductor packages, and all can be installed wtih conda install r-base bioconductor-deseq2 etc.
