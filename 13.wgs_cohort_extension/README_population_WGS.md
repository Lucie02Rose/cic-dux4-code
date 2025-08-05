This folder contains all scripts that were used to investigate 4q DUX4 gene heterogeneity in a population of matched normal patients across projects.
I have downloaded all bam files into a /lustre/scratch126/casm/team274sb/lr26/population_dux/<sample>/<sample>.sample.dupmarked.bam
fashion, so that it is easy to just use sed to replace sample IDs. I have made lists of all samples used in population-xxx-script-generator.sh which are 
scripts used to adapt the central population-dux.sh script and generate script-ID.sh for each sample in the same directory, which I can then submit 
using a for loop, e.g. for script in script_PD*.sh; do bsub < "$script"; done to the farm and script will run whenever cpus are available. 
I have used four types of tumor normals - orphan tumors, sarcoma, Wilms tumor and placenta samples. 
The population-dux.sh script uses minimap2 to subset and extract relevant reads from hg38-aligned bam files and realigns the subset reads to the T2T reference.
The bamscale-dux.sh script then takes all these subset bams and normalises them to FPKM, which is then used in the population_dux_downstream.ipynb script 
to plot and visualise the results.

