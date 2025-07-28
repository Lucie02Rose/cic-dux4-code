#!/bin/bash
#BSUB -n 64
#BSUB -M 400000
#BSUB -R 'span[hosts=1] select[mem>400000] rusage[mem=400000]'
#BSUB -q hugemem
#BSUB -J hifitum
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-tumordenovo.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-tumordenovo.err

### Activate conda environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate hifiasm-env

# List of input FASTQ.gz files
combined_fastq="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_1B01_hifi_reads.fastq.gz"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-tumor"
tmp_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-tumor/tmp"

#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 64
#BSUB -M 400000
#BSUB -R 'span[hosts=1] select[mem>400000] rusage[mem=400000]'
#BSUB -q hugemem
#BSUB -J sequel_hifiasm
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-sequeldenovo.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-sequeldenovo.err

### Acactivate the hifiasm conda environment ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate hifiasm-env

### use the already 30x revio sample ###
combined_fastq="/lustre/scratch126/cellgen/behjati/lr26/PacBio-fastq/tumor_1B01_hifi_reads.fastq.gz"
output_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-revio"
tmp_dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-revio/tmp"

### create the output and temporary directory, export the temporary directory ####
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"
export TMPDIR="$tmp_dir"

### concatenate the sequel II runs to have one 30x depth of coverage file ###
cat "$input_fastq1" "$input_fastq2" "$input_fastq3" > "$combined_fastq"
### run hifiasm on the sequel II ###
hifiasm -o "$output_dir" -t64 "$combined_fastq"
### note that assembly is running ###
echo "assembly process has started"
