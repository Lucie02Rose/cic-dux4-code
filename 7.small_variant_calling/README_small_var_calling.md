This directory contains all scripts related to small variant calling (e.g. SNVs and indels) 
in long read data using PEPPER-Margin-DeepVariant (sif image included in this directory). 
Firstly, variants were called on all three Revio samples (mom, blood, tumor) using both genome references
and all deepvariant- scripts. Then, deepvariant-intersection-somatic.sh was used to filter for somatic variants
in the tumor. After, variant annotation was performed using COSMIC on the hg38 reference-called somatic tumor 
cosmic-hg38.sh, and Variant Effect Predictor (dbSNP and ClinVar) were used to annotate both hg38-aligned 
and t2t-aligned somatic (vep-somatic) and germline (vep-deepvar-blood) variants in the tumor. 
Maternally inherited variants were found by using VEP annotation on the maternal sample. There is no
way of distinguishing paternal variants from de novo variants. There are two separate Python scripts for 
downstream filtering - for cosmic and for VEP. Both mom and tumor VEP filtered outputs are included as Excel files. 
