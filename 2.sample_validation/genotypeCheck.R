### This is a R-based allele counting tool developed by Dr Mi Trinh ###

# Check genotype consistency
# Are all the BAMs you're going to use from the same individual?  Check before you start
# Run Genotype Check for all scRNA-seq data against all WGS data

# Install relevant packages
#devtools::install_github('constantAmateur/alleleIntegrator')
#install.packages('/nfs/users/nfs_m/my4/alleleIntegrator_0.9.1.tar.gz',repos = NULL,type='source')
#BiocManager::install("VariantAnnotation")
#BiocManager::install('SNPRelate')

#sudo apt-get install bcftools
#sudo cp /nfs/users/nfs_m/my4/bin/alleleCounter /usr/local/bin/alleleCounter
#sudo chmod +x /usr/local/bin/alleleCounter

##-------------------##
##   Libraries     ####
##-------------------##
library(tidyverse)
library(alleleIntegrator)
#source('genotypeCheck_helperFunctions.R')

#------------------------ This section must be specified by users -----------------------------------------##

## Set (and create if not exist) output directory
outDir = '/lustre/scratch126/casm/team274sb/lr26/allele-integrator-pbfix'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
setwd(outDir)


#### Get list of RNA and DNA BAM files ####
#----- 10X RNA BAMs
bams10X = c('Ewings-patient' = "/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_bams/GOSH028/NB8113359_possorted_genome_bam.bam",
            'patient tumor FO CG_SB_NB13960948' = "/lustre/scratch126/casm/team274sb/lr26/scRNA/cellranger800_count_47774_CG_SB_NB13960948_GRCh38-1_2_0/possorted_genome_bam.bam",
            'patient tumor FO CG_SB_NB13960949' = "/lustre/scratch126/casm/team274sb/lr26/scRNA/cellranger800_count_47774_CG_SB_NB13960949_GRCh38-1_2_0/possorted_genome_bam.bam",
            'patient tumor FO CG_SB_NB13960950' = "/lustre/scratch126/casm/team274sb/lr26/scRNA/cellranger800_count_47774_CG_SB_NB13960950_GRCh38-1_2_0/possorted_genome_bam.bam",
            'patient tumor FO CG_SB_NB13960951' = "/lustre/scratch126/casm/team274sb/lr26/scRNA/cellranger800_count_47774_CG_SB_NB13960951_GRCh38-1_2_0/possorted_genome_bam.bam",
            'patient tumor FT CG_SB_NB14449539' = "/lustre/scratch126/casm/team274sb/lr26/scRNA/cellranger800_count_48097_CG_SB_NB14449539_GRCh38-1_2_0/possorted_genome_bam.bam",
            'patient tumor FT CG_SB_NB14449540' = "/lustre/scratch126/casm/team274sb/lr26/scRNA/cellranger800_count_48097_CG_SB_NB14449540_GRCh38-1_2_0/possorted_genome_bam.bam",
            'patient tumor FT CG_SB_NB14449541' = "/lustre/scratch126/casm/team274sb/lr26/scRNA/cellranger800_count_48097_CG_SB_NB14449541_GRCh38-1_2_0/possorted_genome_bam.bam")

if(length(bams10X) > 0){
  bams10X = bams10X[file.exists(bams10X)]  
}

# Check that each bam10X file has a unique name
if(length(unique(names(bams10X))) != length(bams10X)){
  stop(sprintf('Duplicated bams10X names detected: %s',names(bams10X)[duplicated(names(bams10X))]))
}

#----- RNA BAMs
rnaBAMs = c('patient-bulk-RNA' = "/lustre/scratch126/casm/team274sb/na15/PacBio/rnaseq/Aligned.sortedByCoord.out.bam")
#names(rnaBAMs) = c()
if(length(rnaBAMs) > 0){
  rnaBAMs = rnaBAMs[file.exists(rnaBAMs)]  
}
# Check that each bam10X file has a unique name
if(length(unique(names(rnaBAMs))) != length(rnaBAMs)){
  stop(sprintf('Duplicated rnaBAMs names detected: %s',names(rnaBAMs)[duplicated(names(rnaBAMs))]))
}

#----- DNA BAMs 
dnaBAMs = c('patient-bulk-DNA' = "/nfs/cancer_ref01/nst_links/live/3306/PD54858d/PD54858d.sample.dupmarked.bam",
            'patient-bulk-DNA-un' = "/nfs/cancer_ref01/nst_links/live/3306/PD54858b/PD54858b.sample.dupmarked.bam",
            'patient-tum-1A01' = "/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274sb-hg38/1_A01/m64094e_230126_154129.hifi_reads_pbmm2-farm22-bam.softclipped.bam",
            'patient-tum-1A02' = "/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274sb-hg38/1_A02/m64178e_230206_134948.hifi_reads_pbmm2-farm22-bam.softclipped.bam",
            'patient-tum-2B01' = "/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274sb-hg38/2_B01/m64178e_230207_165902.hifi_reads_pbmm2-farm22-bam.softclipped.bam",
            'patient-tum-1B01' = "/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274-hg38-revio/1_B01/m84047_230404_172053_s2.hifi_reads.default_pbmm2-farm22-bam.softclipped.bam",
            'patient-blood-1B02' = "/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274-hg38-revio/1_B02/m84047_240202_152510_s2.hifi_reads.bc2025_pbmm2-farm22-bam.softclipped.bam",
            'mum-blood-1C01' = "/lustre/scratch126/casm/team274sb/lr26/pbmm2-alignment-team274-hg38-revio/1_C01/m84047_240202_155616_s3.hifi_reads.bc2026_pbmm2-farm22-bam.softclipped.bam",
            'mum-bulk-DNA' = "/nfs/cancer_ref01/nst_links/live/3306/PD54859b/PD54859b.sample.dupmarked.bam" 
           )
#names(dnaBAMs) = c()
if(length(dnaBAMs) > 0){
  dnaBAMs = dnaBAMs[file.exists(dnaBAMs)]
}

# Check that each dnaBAM file has a unique name
if(length(unique(names(dnaBAMs))) != length(dnaBAMs)){
  stop(sprintf('Duplicated dnaBAMs names detected: %s',names(dnaBAMs)[duplicated(names(dnaBAMs))]))
}
#------------------------ End of essential user-specified section -----------------------------------------##



#--------------------------------------------------------- #
# Parameters - users can adapt these parameters if needed #
#--------------------------------------------------------- #
wgs_ref_version = 'hg38' # or 'hg19'
sc_ref_version = 'hg38' # or 'hg19'
refGenome=NULL
refGenome10X=NULL

# reference genome used for WGS
if(is.null(refGenome)){
  if(wgs_ref_version == 'hg38'){
    refGenome = '/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa' # hg38  
  }else{
    refGenome = '/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh37d5/genome.fa' #hg19
  }
}

if(is.null(refGenome10X)){
  if(sc_ref_version == 'hg38'){
    # reference genome used for scRNA-seq
    refGenome10X = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/fasta/genome.fa' # hg38  
    gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf' # hg38
  }else{
    # This is casm hg19 reference genome
    refGenome10X = '/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCH37d5/genome.fa'  # hg19
    gtf = '/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh37d5/star/e75/ensembl.gtf' # hg19
  }
}


liftChain = '/lustre/scratch125/casm/team274sb/mt22/hg19ToHg38_noChr.over.chain'
verbose = TRUE
nParallel = 48
downsampleBams = FALSE
downsampleBamFraction=0.01
removeBam = FALSE
skipIfExists = TRUE
nMaxSim = 8 # Max number of samples being processed at once
plotDir = outDir

if(!dir.exists(file.path(outDir, "GenotypingResults"))){
  dir.create(file.path(outDir, "GenotypingResults"), recursive = TRUE)
}

# Create output directory
if(!dir.exists(outDir)){
  dir.create(file.path(outDir,'downsampledBams'),recursive = T)
  dir.create(file.path(outDir,'GenotypingResults'),recursive = T)
}


#############################

if(downsampleBams){
  ds_bams = downsampleBam(bams = c(dnaBAMs,rnaBAMs,bams10X),
                          outDir = file.path(outDir,'downsampledBams'),
                          genotyping_outDir = file.path(outDir,'GenotypingResults'),
                          downsampleBamFraction=downsampleBamFraction,
                          nParallel=nParallel,skipIfExists=skipIfExists)
  
  genotyping_outputs = file.path(outDir,'GenotypingResults',paste0(names(ds_bams),'_genotypeCheck.tsv'))
  BAMs = ds_bams
}else{
  genotyping_outputs = file.path(outDir,'GenotypingResults',paste0(c(names(dnaBAMs),names(rnaBAMs),names(bams10X)),'_genotypeCheck.tsv'))
  BAMs = c(dnaBAMs,rnaBAMs,bams10X)
}



genoCheck = matchBAMs(BAMs = BAMs,
                      refGenomes = rep(c(refGenome,refGenome,refGenome10X),c(length(dnaBAMs),length(rnaBAMs),length(bams10X))),
                      outputs = genotyping_outputs,
                      liftOvers=rep(c(liftChain,liftChain,liftChain),c(length(dnaBAMs),length(rnaBAMs),length(bams10X))),
                      is10X=rep(c(FALSE,FALSE,TRUE),c(length(dnaBAMs),length(rnaBAMs),length(bams10X))),
                      nParallel=nParallel,nMaxSim=nMaxSim,nChunks=6,skipIfExists=skipIfExists)

if(!is.null(plotDir)){
  library(RColorBrewer)
  library(circlize)
  library(ComplexHeatmap)
  bamGrouping = NULL
  colPal = 'Greens'
  colFun = suppressWarnings(brewer.pal(100, colPal))
  colFun = colorRamp2(seq(0.5, 1, length.out = length(colFun)), 
                      colFun)
  hm = Heatmap(genoCheck$ibs$ibs, col = colFun, name = "IBS", show_row_names = TRUE, 
               show_column_names = TRUE, show_row_dend = FALSE, 
               show_column_dend = FALSE, row_title_rot = 0, column_split = bamGrouping, 
               row_split = bamGrouping)
  pdf(paste0(plotDir,'/genotypeCheck_',downsampleBamFraction,'dsDepth.pdf'))
  draw(hm)
  dev.off()
  
}
#If anything is less than 0.8 and you should be concerned...
message(sprintf("The minimum similarity found was %g",min(genoCheck$ibs$ibs)))


if(removeBam){
  message('Deleting downsampled BAM files')
  system(sprintf("rm -r %s", file.path(outDir,'downsampledBams',names(c(dnaBAMs,rnaBAMs,bams10X)),'_',downsampleBamFraction,'_downsampled.bam')))
}
