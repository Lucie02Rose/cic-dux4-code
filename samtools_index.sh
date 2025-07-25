module load samtools

cd /lustre/scratch126/casm/team274sb/na15/PacBio/ 
samtools sort CICSARC.merged.bam -o CICSARC.merged_sort.bam
samtools index CICSARC.merged_sort.bam
