This folder contains scripts connected to quality control, alignment of long-read data
to the T2T and hg38 references, indexing and optional decompression of the references.

All scripts work with the chm13v2.0.fa reference and assume that it has been downloaded already using 
the wget commands listed in the README_references.md in references. 
chm13v2.0.fa can either be downloaded compressed or decompressed (there is an optional decompressing step but I downloaded the already decompressed file). Note that the indexing.sh scripts need to be run first for the alignments to work.

There are files for individual samples, or if there is more compute available, for all samples in a for loop.
The bam_to_fastq file contains conversion of samples to fastq, which is important for combining all the tumor
runs together (tumor_alignment_combined) and also for the de novo assemblies highlighted in the de_novo_assembly section. 

The base conda environment is already included in my .bashrc, so I technically do not need to activate
it. However, the content of the conda environment matters (note that it uses the free version
of anaconda, e.g. was installed with miniforge and not miniconda) since anaconda is now paid. 

The programms in the base environment can be found in the .yaml file also in this folder.
The pbmm2 version used is also listed in the .yaml. Note that pbmm2 is indifferent if the index
is named .fa.mmi or only .mmi.
