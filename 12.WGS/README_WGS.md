This directory contains WGS outputs. Firstly, there is a wgs-alignment-t2t.sh script which 
is used to extract fastq reads from bam files and align the reads to the T2T reference. 
Then, there is a WGS.ipynb script used for downstread filtering of variant caller outputs 
from the CASM-IT pipeline, of which the filtered outputs are included in this directory as well. 
There is an additional symlinks.R script that can be used internally to search for sample outputs
throughout the file system. 
