#BSUB -n 16
#BSUB -M 200000
#BSUB -R 'span[hosts=1] select[mem>200000] rusage[mem=200000]'
#BSUB -q basement
#BSUB -J ngmlr-alignment-hg38-revio-1B01
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### Directories
reference="/lustre/scratch126/casm/team274sb/lr26/hg38/hg38.fna"
input_fastq="/lustre/scratch126/casm/team274sb/lr26/fastq-revio/m84047_230404_172053_s2.hifi_reads.default.fastq.gz"
output_bam="/lustre/scratch126/casm/team274sb/lr26/ngmlr-alignment-hg38/1_B01-revio/m84047_230404_172053_s2.hifi_reads.default.bam"

### Run NGMLR alignment

echo "Aligning $input_bam to $reference..."
ngmlr -r "$reference" -q "$input_fastq" -o "$output_bam" -t 16 -x pacbio

echo "NGMLR alignment completed: $output_bam"

