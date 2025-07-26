#!/bin/bash
#BSUB -n 16
#BSUB -M 120000
#BSUB -R 'span[hosts=1] select[mem>120000] rusage[mem=120000]'
#BSUB -q basement
#BSUB -J liftover
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/liftover-l.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/liftover-l.err

### Activate Conda Environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
cosmic="/lustre/scratch125/casm/team274sb/na15/cosmic/Cosmic_MutantCensus_v99_GRCh38.tsv.gz"
chain="/lustre/scratch126/casm/team274sb/lr26/T2T/grch38-chm13v2.chain"
T2T="/lustre/scratch126/casm/team274sb/lr26/T2T"
liftover="/software/cellgen/team274/lr26/miniforge3/bin/liftOver"

cd "$T2T"

zcat "$cosmic" | \
awk -F'\t' '
BEGIN {OFS="\t"}
NR==1 {
  # Mapping headers
  for (i=1; i<=NF; i++) header[$i] = i
  next
}
{
  # Extract the necessary columns
  chr = $header["CHROMOSOME"]
  start = $header["GENOME_START"]
  end = $header["GENOME_STOP"]
  id = $header["MUTATION_ID"]

  # Check if necessary columns are not empty and valid (numbers for start/end)
  if (chr != "" && start ~ /^[0-9]+$/ && end ~ /^[0-9]+$/ && id != "") {
    print "chr"chr, start - 1, end, id  # Convert to zero-based for liftOver
  }
}' > cosmic.hg38.liftover.bed


"$liftover" cosmic.hg38.liftover.bed "$chain" cosmic.t2t.liftover.bed cosmic.unmapped.liftover.bed

awk -F'\t' 'BEGIN{OFS="\t"} {print $4, $1, $2+1, $3}' cosmic.t2t.liftover.bed > id_to_t2t_coords.liftover.tsv

zcat "$cosmic" | \
awk -F'\t' 'BEGIN{OFS="\t"}
FNR==NR {
  split($0, fields, "\t")
  coord_map[fields[0]] = fields[1] "\t" fields[2] "\t" fields[3]
  next
}
NR==1 {
  for (i=1; i<=NF; i++) {
    header[i]=$i
    if ($i == "MUTATION_ID") id_col = i
    if ($i == "CHROMOSOME") chr_col = i
    if ($i == "GENOME_START") start_col = i
    if ($i == "GENOME_STOP") stop_col = i
  }
  print
  next
}
{
  id = $id_col
  if (id in coord_map) {
    split(coord_map[id], coords, "\t")
    $chr_col = coords[0+1]
    $start_col = coords[1+1]
    $stop_col = coords[2+1]
  }
  print
}' id_to_t2t_coords.liftover.tsv - > Cosmic_MutantCensus_v99_T2T.liftover.tsv














