

# Step 1: Intersect VCFs to retain tumor-specific variants
bcftools isec -c none -n -1 -w1 -O v \
-o tumor_only.vcf.gz \
~/my_lustre/sawfish-tumor-all/tumor-all_pbmm2-farm22_joint_call/genotyped.sv.vcf.gz \
~/my_lustre/sawfish/1B02genotyped.sv.vcf.gz

# Index the output
tabix -p vcf tumor_only.vcf.gz

# Step 2: Filter out variants with >5% allele frequency in normal
bcftools filter -e 'AD[1:1] / (AD[1:0] + AD[1:1]) > 0.05' tumor_only.vcf.gz -O v -o tumor_somatic.vcf.gz

# Index the filtered file
tabix -p vcf tumor_somatic.vcf.gz

# Step 3: Keep only large structural variants (SVLEN ≥ 1000)
bcftools filter -i 'INFO/SVLEN>=1000' tumor_somatic.vcf.gz -O v -o tumor_large_svs.vcf.gz

# Index the final file
tabix -p vcf tumor_large_svs.vcf.gz

# Step 4: Extract data for visualization
bcftools query -f '%CHROM\t%POS\t%SVTYPE\t%MATEID\n' tumor_large_svs.vcf.gz > sv_for_plot.tsv
