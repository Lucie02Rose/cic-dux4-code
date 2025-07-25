module load samtools

cd /lustre/scratch126/casm/team274sb/na15/PacBio/NGMLR/

/lustre/scratch126/casm/team274sb/na15/SComatic/python/miniconda/envs/SComatic/bin/sniffles --input /lustre/scratch126/casm/team274sb/na15/PacBio/alignment/CICSARC.merged_sort.bam -v CICSARC.merged_pbmm2_wtr.vcf --non-germline
