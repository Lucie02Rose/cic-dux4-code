#!/bin/bash
#BSUB -n 16
#BSUB -M 120000
#BSUB -R 'span[hosts=1] select[mem>120000] rusage[mem=120000]'
#BSUB -q basement
#BSUB -J liftover
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/liftover.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/liftover.err

### Activate Conda Environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
cosmic="/lustre/scratch125/casm/team274sb/na15/cosmic/Cosmic_MutantCensus_v99_GRCh38.tsv.gz"
chain="/lustre/scratch126/casm/team274sb/lr26/T2T/grch38-chm13v2.chain"
T2T="/lustre/scratch126/casm/team274sb/lr26/T2T"

cd "$T2T"

zcat "$cosmic" | \
awk -F'\t' '
BEGIN {OFS="\t"}
NR==1 {
  for (i=1; i<=NF; i++) col[$i]=i
  next
}
{
  chr="chr"$col["CHROMOSOME"];
  start=$col["GENOME_START"]-1;  # BED is 0-based
  end=$col["GENOME_STOP"];
  id=$col["MUTATION_ID"];
  print chr, start, end, id
}' > cosmic.hg38.bed

CrossMap bed "$chain" cosmic.hg38.bed "$reference" cosmic.t2t.bed

awk -F'\t' 'BEGIN{OFS="\t"} {print $4, $1, $2+1, $3}' cosmic.t2t.bed > id_to_t2t_coords.tsv

zcat "$cosmic" | \
awk -F'\t' 'BEGIN{OFS="\t"}
FNR==NR {
  coord_map[$1]=$2"\t"$3"\t"$4;  # MUTATION_ID → chr \t start \t end
  next
}
NR==1 {
  for (i=1; i<=NF; i++) header[i]=$i
  for (i=1; i<=NF; i++) {
    if (header[i] == "MUTATION_ID") id_col = i
    if (header[i] == "CHROMOSOME") chr_col = i
    if (header[i] == "GENOME_START") start_col = i
    if (header[i] == "GENOME_STOP") stop_col = i
  }
  print
}
NR>1 {
  id = $id_col
  if (id in coord_map) {
    split(coord_map[id], coords, "\t")
    $chr_col = coords[1]
    $start_col = coords[2]
    $stop_col = coords[3]
  }
  print
}' id_to_t2t_coords.tsv - > Cosmic_MutantCensus_v99_T2T.tsv














