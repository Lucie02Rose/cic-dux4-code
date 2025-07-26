This folder contains scipts connected to alignment of long-read data
to the T2T reference, indexing and optional decompression of the reference.

All scripts work with the chm13v2.0.fa reference and assume that it has been downloaded already using 
the wget commands listed in the README_references.md. chm13v2.0.fa can either be downloaded compressed
or decompressed (there is an optional decompressing step but I downloaded the already decompressed 
file). Note that the indexing.sh script needs to be run first for the alignments to work.
Alignment files are organised by runs - e.g. tumor, blood, mom.

The base conda environment is already included in my .bashrc, so I technically do not need to activate
it. However, the content of the conda environment matters (note that it uses the free version
of anaconda, e.g. was installed with miniforge and not miniconda) since anaconda is now paid. 

The programms in the base environment can be found in the .yml file also in this folder.
The pbmm2 version used is also listed in the .yml. Note that pbmm2 is indifferent if the index
is named .fa.mmi or only .mmi.
