#!/bin/bash
#BSUB -n 64
#BSUB -M 32000
#BSUB -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q long
#BSUB -J cnv-pacbio
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-cnv-pacbio.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-cnv-pacbio.err

### Activate Conda Environment
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### Directories
ref="/lustre/scratch126/cellgen/behjati/lr26/T2T/"
reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
dbsnp="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_dbSNPv155.vcf.gz"
tumor="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/tumor_all_4_hifi_reads_pbmm2.bam"
mom="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/mom_1B02_hifi_reads_pbmm2.bam"
patient="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/1C01_blood_hifi_reads_pbmm2.bam"
dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-hificnv-tumor-all"
annot="/lustre/scratch126/cellgen/behjati/lr26/T2T/sorted_positions.tsv"
cd "$ref"

#sample only dbsnps on chr1 and 10 and index the vcf
# 1. Extract all positions (chrom and pos) from the entire dbSNP (no filtering here)

#bcftools view -H "$dbsnp" | cut -f1,2 > all_positions.tsv
# 2. Sort + unique to get unique positions genome-wide
#sort all_positions.tsv | uniq > unique_positions.tsv
# 3. Keep every 100th position genome-wide
#awk 'NR % 100 == 0' unique_positions.tsv > subsampled_positions.tsv
# 4. Filter subsampled positions to keep only chr1 and chr10
#grep -P '^(chr1|chr10)\t' subsampled_positions.tsv > positions_100th_unique_chr1_10.tsv

sort -k1,1V -k2,2n subsampled_positions.tsv > "$annot" 

cd "$dir"

split -l 50000 "$annot" chunk_

run_mpileup() {
    local sample_prefix=$1
    local bamfile=$2

    echo "Processing $sample_prefix..."

    # Use a unique prefix for chunk outputs per sample to avoid collisions
    parallel -j 64 bcftools mpileup -f "$reference" -a FORMAT/AD,FORMAT/DP -R {} -o "${sample_prefix}_{}.bcf" "$bamfile" ::: chunk_*

    # Merge all chunk BCFs for this sample
    bcftools concat -O b -o merged_${sample_prefix}.bcf ${sample_prefix}_*.bcf

    # Remove chunk BCFs for this sample to clean up
    rm -f ${sample_prefix}_*.bcf

    # Convert merged BCF to VCF
    bcftools view merged_${sample_prefix}.bcf -Ov -o merged_${sample_prefix}.vcf

    # Extract and calculate allele fractions
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t[%DP]\n' merged_${sample_prefix}.vcf | \
    awk -F'\t' '{
        split($5, ad, ",");
        total=0;
        for (i in ad) total+=ad[i];
        if (total > 0) {
            printf "%s", $0;
            for (i in ad) printf "\t%.3f", ad[i] / total;
            print "";
        } else {
            print $0, "0";
        }
    }' > "${sample_prefix}_fin_new_all_merged.vcf"

    echo "$sample_prefix done."
}

# Run for tumor
run_mpileup patient "$tumor"

rm chunk_*
