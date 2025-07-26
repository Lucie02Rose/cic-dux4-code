# this script finds files in irods and makes symlinks to all files 
# of a PD sample given the target directory, project number and PD ID
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Validate input - activate environment with conda activate my-r an then run R script
# have to input target directory and project number and patient ID
if (length(args) < 3) {
  stop("Usage: Rscript symlinks.R <target_dir> <project_number> <pd_id>")
}

target_dir <- args[1]
project_number <- args[2]  
pd_id <- args[3]           

# Create target_dir if it doesn't exist
if (!dir.exists(target_dir)) {
  dir.create(target_dir, recursive = TRUE)
}

setwd("/nfs/")

# List all directories matching "irods-cgp"
irods <- list.dirs(path = ".", recursive = FALSE, full.names = TRUE)
irods <- grep("irods-cgp", irods, value = TRUE)

# Find BAM files and create symlinks
nfs_paths_list <- lapply(irods, function(irods_path) {
  sample_path <- file.path(irods_path, paste0("intproj/", project_number, "/sample"))

  if (!dir.exists(sample_path)) return(character(0))

  PDs <- list.dirs(path = sample_path, recursive = FALSE, full.names = TRUE)
  PDs <- grep(pd_id, PDs, value = TRUE)

  bam_files <- unlist(lapply(PDs, function(dir) {
    files <- list.files(path = dir, pattern = pd_id, recursive = TRUE, full.names = FALSE)
    file.path(normalizePath(dir), files)
  }))

  for (bam_file in bam_files) {
    link_name <- file.path(target_dir, basename(bam_file))

    if (!file.exists(link_name)) {
      success <- try(file.symlink(bam_file, link_name), silent = TRUE)
      if (inherits(success, "try-error")) {
        warning(paste("Failed to create symlink for:", bam_file))
      }
    }
  }

  bam_files
})

nfs_paths <- unlist(nfs_paths_list, use.names = FALSE)
print(nfs_paths) 
# path lists can also be obtained by ls -lh in the target directory after linking
