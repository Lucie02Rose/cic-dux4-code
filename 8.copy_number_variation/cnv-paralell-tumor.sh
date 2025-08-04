#!/bin/bash
### parameters for LSF job ###
#BSUB -n 64
#BSUB -M 32000
#BSUB -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]'
#BSUB -q long
#BSUB -J cnv-pacbio
#BSUB -G team274
#BSUB -o /lustre/scratch126/cellgen/behjati/lr26/outputs/%J-cnv-pacbio-tumor.out
#BSUB -e /lustre/scratch126/cellgen/behjati/lr26/errors/%J-cnv-pacbio-tumor.err

### activate the environment with latest bcftools ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate bioinfo

### define all directories and files ###
ref="/lustre/scratch126/cellgen/behjati/lr26/T2T/"
reference="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0.fa"
### here I am using the dbSNP database lifted over to the T2T reference ###
dbsnp="/lustre/scratch126/cellgen/behjati/lr26/T2T/chm13v2.0_dbSNPv155.vcf.gz"
tumor="/lustre/scratch126/cellgen/behjati/lr26/PacBio-aligned/tumor_all_4_hifi_reads_pbmm2.bam"
dir="/lustre/scratch126/cellgen/behjati/lr26/PacBio-hificnv"
annot="/lustre/scratch126/cellgen/behjati/lr26/T2T/sorted_positions.tsv"

### change to the reference directory ###
cd "$ref"

### the dbSNP vcf needs to be sampled for positions ###
### extract all positions (just the first two fields which are chrom and pos) from the entire dbSNP ###
bcftools view -H "$dbsnp" | cut -f1,2 > all_positions.tsv
### then sort the file and get all the unique positions genome-wide ###
### (e.g. there are several SNPs possible for 1 position - ex. A, C, G, T) but here we don't care ###
### and don't want to sample the same locus multiple times ###
sort all_positions.tsv | uniq > unique_positions.tsv
### to make it faster, subset to every 100th unique position ###
awk 'NR % 100 == 0' unique_positions.tsv > subsampled_positions.tsv
### another option is to then filter for chromosomes of interest - not done here ###
#grep -P '^(chr1|chr10)\t' subsampled_positions.tsv > positions_100th_unique_chr1_10.tsv
### sort the subsampled file accordingly ###
sort -k1,1V -k2,2n subsampled_positions.tsv > "$annot" 

### now change to the output directory or make it if it is not yet made ###
mkdir -p "$dir"
cd "$dir"

### split the annotation file into smaller files each of 50 000 lines ###
split -l 50000 "$annot" chunk_

### this function runs the bcftools mpileup using the sample prefix and the bam file ###
run_mpileup() {
    local sample_prefix=$1
    local bamfile=$2
    ### prints which sample is processed ###
    echo "Processing $sample_prefix..."
    ### uses parallel with 64 cpus and a unique prefix for each annot chunk output per sample to avoid collisions ###
    ### uses also the T2T reference and we are interested in the allele depth and sequencing depth ###
    ### using bcf since it is lighter ###
    parallel -j 64 bcftools mpileup -f "$reference" -a FORMAT/AD,FORMAT/DP -R {} -o "${sample_prefix}_{}.bcf" "$bamfile" ::: chunk_*
    ### this concatenates all the bcf chunks into one final merged one ###
    bcftools concat -O b -o merged_${sample_prefix}.bcf ${sample_prefix}_*.bcf
    ### now get rid of the chunk bcfs because there are not needed and take up the file quota ###
    rm -f ${sample_prefix}_*.bcf
    ### convert bcf to the standard vcf ###
    bcftools view merged_${sample_prefix}.bcf -Ov -o merged_${sample_prefix}.vcf
    ### this extracts and calculates the allele fractions depending on depth at that locus and adds to the final vcf ###
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
    ### prints out that it is done ###
    echo "$sample_prefix done."
}

### use this function on the sample defined above (e.g. tumor)
run_mpileup tumor "$tumor"
### remove the chunk files ###
rm chunk_*
