#!/bin/bash
### parameters for the LSF job ###
#BSUB -n 16
#BSUB -M 100000
#BSUB -R 'span[hosts=1] select[mem>100000] rusage[mem=100000]'
#BSUB -q basement
#BSUB -J improv_allele-bcftools
#BSUB -G team274
#BSUB -o /lustre/scratch126/casm/team274sb/lr26/output_logs/%J.out
#BSUB -e /lustre/scratch126/casm/team274sb/lr26/error_logs/%J.err

### conda environment activation (new versions of bcftools and samtools in base) ###
source /software/cellgen/team274/lr26/miniforge3/bin/activate
conda activate base

### define directories for the reference, positions to be sampled, all samples ###
reference="/lustre/scratch126/casm/team274sb/lr26/hg38/genome.fa"
positions_bed="/lustre/scratch126/casm/team274sb/lr26/allele-integrator-pbfix/GenotypingResults/patient_output.bed"
s1A01="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274sb-hg38/1_A01/m64094e_230126_154129.hifi_reads_pbmm2-farm22-bam.bam"
s1A02="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274sb-hg38/1_A02/m64178e_230206_134948.hifi_reads_pbmm2-farm22-bam.bam"
s2B01="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274sb-hg38/2_B01/m64178e_230207_165902.hifi_reads_pbmm2-farm22-bam.bam"
s1B01="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274-hg38-revio/1_B01/m84047_230404_172053_s2.hifi_reads.default_pbmm2-farm22-bam.bam"
s1B02="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274-hg38-revio/1_B02/m84047_240202_152510_s2.hifi_reads.bc2025_pbmm2-farm22-bam.bam"
s1C01="/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274-hg38-revio/1_C01/m84047_240202_155616_s3.hifi_reads.bc2026_pbmm2-farm22-bam.bam"
mom="/nfs/cancer_ref01/nst_links/live/3306/PD54859b/PD54859b.sample.dupmarked.bam"
patient="/nfs/cancer_ref01/nst_links/live/3306/PD54858d/PD54858d.sample.dupmarked.bam"

### for each of the samples defined above, calculate the number of nucleotides at the bed file position compared to reference
### save allele depth information depending on the no of allele counts at the location, store in a bcf file, convert to vcf 
### use awk to calculate total depth at each position and calculate alelle proportion dividing each by total
### output the original and add the calculated allele frequency proportions to the new vcf

### do this for all samples ###
bcftools mpileup -f "$reference" -a FORMAT/AD -R "$positions_bed" -o s1A01.bcf "$s1A01"
bcftools view s1A01.bcf -Ov -o s1A01.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' s1A01.vcf | \
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
}' > s1A01_fin.vcf

bcftools mpileup -f "$reference" -a FORMAT/AD -R "$positions_bed" -o s1A02.bcf "$s1A02"
bcftools view s1A02.bcf -Ov -o s1A02.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' s1A02.vcf | \
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
}' > s1A02_fin.vcf

bcftools mpileup -f "$reference" -a FORMAT/AD -R "$positions_bed" -o s2B01.bcf "$s2B01"
bcftools view s2B01.bcf -Ov -o s2B01.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' s2B01.vcf | \
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
}' > s2B01_fin.vcf

bcftools mpileup -f "$reference" -a FORMAT/AD -R "$positions_bed" -o s1B01.bcf "$s1B01"
bcftools view s1B01.bcf -Ov -o s1B01.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' s1B01.vcf | \
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
}' > s1B01_fin.vcf

bcftools mpileup -f "$reference" -a FORMAT/AD -R "$positions_bed" -o s1B02.bcf "$s1B02"
bcftools view s1B02.bcf -Ov -o s1B02.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' s1B02.vcf | \
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
}' > s1B02_fin.vcf

bcftools mpileup -f "$reference" -a FORMAT/AD -R "$positions_bed" -o s1C01.bcf "$s1C01"
bcftools view s1C01.bcf -Ov -o s1C01.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' s1C01.vcf | \
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
}' > s1C01_fin.vcf

bcftools mpileup -f "$reference" -a FORMAT/AD -R "$positions_bed" -o patient.bcf "$patient"
bcftools view patient.bcf -Ov -o patient.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' patient.vcf | \
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
}' > patient_fin.vcf

bcftools mpileup -f "$reference" -a FORMAT/AD -R "$positions_bed" -o mom.bcf "$mom"
bcftools view mom.bcf -Ov -o mom.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' mom.vcf | \
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
}' > mom_fin.vcf















