module load samtools

cd /lustre/scratch126/casm/team274sb/na15/PacBio/alignment/
samtools flagstat CICSARC.merged_sort.bam > CICSARC.merged_sort.flagstat
samtools flagstat m64094e_230126_154129_sort.bam > m64094e_230126_154129_sort.flagstat
samtools flagstat m64178e_230206_134948_sort.bam > m64178e_230206_134948_sort.flagstat
samtools flagstat m64178e_230207_165902_sort.bam > m64178e_230207_165902_sort.flagstat
