This directory contains all GIAB scripts, from alignment of all raw bam files
using pbmm2, then conversion of raw files to fastq and their concatenation (giab-to-fastq-and-concat.sh)
which is necessary prequisite for all HiFiASM assemblies. After assemblies, the minimap-giab-hifiasm.sh
script is used to convert contigs to fasta and then align back to the reference. 
Fasta ouputs are also used in the blast-GIAB.sh to carry out a BLAST of four regions of interest in the
4q DUX4 area. Note that the fasta files of teh four regions of interest are in the 9.de_novo_assembly directory.
