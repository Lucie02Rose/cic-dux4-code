module load samtools

cd /lustre/scratch126/casm/team274sb/na15/PacBio/shells/ 
samtools merge -r -b /lustre/scratch126/casm/team274sb/na15/PacBio/shells/bamlist -o /lustre/scratch126/casm/team274sb/na15/PacBio/CICSARC.merged.bam
