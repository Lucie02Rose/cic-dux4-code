module load samtools

cd /lustre/scratch126/casm/team274sb/na15/PacBio/NGMLR/

/lustre/scratch126/casm/team274sb/na15/SComatic/python/miniconda/envs/SComatic/bin/sniffles --input /lustre/scratch126/casm/team274sb/na15/PacBio/alignment/m64094e_230126_154129_sort.bam -v m64094e_230126_154129_pbmm2_wtr.vcf --non-germline
