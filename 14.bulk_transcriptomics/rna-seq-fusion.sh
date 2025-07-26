#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.o
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.e
#BSUB -n 16
#BSUB -M 300000
#BSUB -R 'span[hosts=1] select[mem>300000] rusage[mem=300000]'
#BSUB -q hugemem
#BSUB -J rna-t2t-alignemnt
#BSUB -G team274

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate star_env

### Directories to process
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
input_dir="/lustre/scratch126/casm/team274sb/lr26/bulk-rna-t2t"

cd "$input_dir"

STAR-Fusion \
   --genome_lib_dir "$reference"
   --left_fq reads_R1.fastq \
   --right_fq reads_R2.fastq \
   --output_dir fusion \
