#!/bin/bash
#BSUB -n 16
#BSUB -M 120000
#BSUB -R 'span[hosts=1] select[mem>120000] rusage[mem=120000]'
#BSUB -q basement
#BSUB -J liftover
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/cosmic.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/cosmic.err

### Activate Conda Environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### Directories
reference="/lustre/scratch126/casm/team274sb/lr26/T2T/chm13v2.0.fa"
cosmic="/lustre/scratch126/casm/team274sb/lr26/T2T/Cosmic_MutantCensus_v99_T2T.liftover.fixed.tsv"
output="/lustre/scratch126/casm/team274sb/lr26/T2T/Cosmic_MutantCensus_v99_T2T.liftover.vcf"
T2T="/lustre/scratch126/casm/team274sb/lr26/T2T"

cd "$T2T"

# Create the VCF header
echo "##fileformat=VCFv4.2" > "$output"
echo "##reference=$reference" >> "$output"
echo "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" >> "$output"
echo "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">" >> "$output"
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO" >> "$output"

# Read the TSV and convert each line to VCF format
awk -F'\t' '
BEGIN {OFS="\t"}
NR==1 {next}  # skip header
{
  chrom = $15  # CHROMOSOME
  pos = $16  # GENOME_START
  id = $7  # COSMIC Variant ID like COSV...
  ref = $24  # GENOMIC_WT_ALLELE
  alt = $25  # GENOMIC_MUT_ALLELE

  gene_symbol = $1
  zygosity = $13
  mut_desc = $11
  mut_cds = $10
  mut_aa = $12
  pubmed = $18
  study_id = $19
  transcript = $3
  sample = $5

  info = "GENE_SYMBOL=" gene_symbol \
  ";MUTATION_ZYGOSITY=" zygosity \
  ";MUTATION_DESCRIPTION=" mut_desc \
  ";MUTATION_CDS=" mut_cds \
  ";MUTATION_AA=" mut_aa \
  ";PUBMED_PMID=" pubmed \
  ";COSMIC_STUDY_ID=" study_id \
  ";TRANSCRIPT_ACCESSION=" transcript \
  ";SAMPLE_NAME=" sample

  print chrom, pos, id, ref, alt, ".", ".", info
}' "$cosmic" >> "$output"












