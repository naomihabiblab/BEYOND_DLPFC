## Setup of CellBender V0.2.0 on machine
## Based on installation instructions from https://github.com/broadinstitute/CellBender, date: 01.03.21
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
# bash ~/miniconda.sh -b -p $HOME/miniconda -f
# 
# ~/miniconda/bin/conda create -n cellbender python=3.7
# source ~/miniconda/bin/activate cellbender
# ~/miniconda/bin/conda install -c anaconda pytables
# ~/miniconda/bin/conda install pytorch torchvision -c pytorch
# pip install -e CellBender


# -------------------------------------------------- #
# Receive command-line arguments for run             #
# -------------------------------------------------- #
suppressPackageStartupMessages(library(optparse))
args <- parse_args(OptionParser(usage="%prog [options]", 
                                option_list=list(
                                  make_option(c("-p","--path"), 
                                              default = "1. Library preprocessing/data/snRNA-seq libraries",
                                              action = "store", 
                                              help="Path to library folders"),
                                  
                                  make_option(c("-b","--bender.folder"), 
                                              default = "CellBender",
                                              action = "store", 
                                              help="Name of folder for storing and retrieving CellBender results"),
                                  
                                  make_option(c("-a","--analysis.folder"), 
                                              default = "Analysis",
                                              action = "store",
                                              help="Name of folder for storing and retrieving analysis results"),
                                  
                                  make_option(c("-j", "--job.details"), 
                                              default = "1. Library preprocessing/data/library.preparation.status.xlsx",
                                              action="store",
                                              help="Path to file containing progress status for each library"),
                                  
                                  make_option(c("-e", "--execute"), 
                                              type = "logical",
                                              action="store",
                                              help="Remove old files and create new ones.")
                                  )))
unloadNamespace("optparse")


# -------------------------------------------------- #
# Update library preparation status and jobs to run  #
# -------------------------------------------------- #
details <- openxlsx::read.xlsx(args$job.details, rowNames = T)
unmanaged.folders <- c(); bender.run <- c(); analysis.run <- c(); missing.demux <- c();
libs <- list.dirs(args$path, full.name=F, recursive = F)

library(rvest)
for(lib in libs) {
  # Examine only libraries managed by the details spreadsheet
  if(! lib %in% rownames(details)) {
    next
  }
  
  # Extract expected number of cells from library's CellRanger report
  if(is.na(details[lib,]$expected.cells) & file.exists(file.path(args$path, lib, 'web_summary.html'))) {
    details[lib,]$expected.cells <- 
      as.numeric(gsub(",", "",jsonlite::fromJSON(gsub("\n\\s+?$", "", gsub(".*const data = ", "", 
           html_text((read_html(file.path(args$path, lib, 'web_summary.html')) %>% html_nodes("script"))[1]))))[["summary"]][["summary_tab"]][["filtered_bcs_transcriptome_union"]][["metric"]]))
  }

  
  # Check existence of output files
  ranger.output.exists   <- file.exists(file.path(args$path, lib, 'raw_feature_bc_matrix.h5'))
  bender.output.exists   <- file.exists(file.path(args$path, lib, args$bender.folder, paste0(lib, "_filtered.h5")))
  demux.output.exists    <- file.exists(file.path(args$path, lib, "demultiplexed_allsites_gt.best"))
  analysis.output.exists <- file.exists(file.path(args$path, lib, args$analysis.folder, paste0(lib, ".seurat.rds")))
  
  
  # Update library's output file status
  details[lib,]$ranger.complete   <- ranger.output.exists
  details[lib,]$bender.complete   <- bender.output.exists
  details[lib,]$demux.complete    <- demux.output.exists
  details[lib,]$analysis.complete <- analysis.output.exists
  
  # If no CellRanger output stop processing library
  if(!ranger.output.exists) next
  
  # If library is supposed to be running:
  #   If still running do nothing
  #   If not running  state Completed/Failed according to existence of output
  if(!is.na(details[lib,]$job.id)) {
    is.running <- !identical(as.character(suppressWarnings(system(paste("squeue | grep", details[lib,]$job.id), intern = T))), character(0))
    if(is.running) {
      # Job is still running
      details[lib,]$job.status <- "Running"
    } else {
      # Logged job.id isn't in queue - completed or failed
      if(details[lib,]$job.name == "CellBender") { 
        details[lib,]$job.status <- ifelse(bender.output.exists, "Completed", "Failed")
        
      } else { if(details[lib,]$job.name == "LibraryAnalysis") {
        details[lib,]$job.status <- ifelse(analysis.output.exists, "Completed", "Failed")
      }}
    }
  }

  # According to running status set jobs to be queued
  if(is.na(details[lib,]$job.status)) {
    # If Ranger and Bender outputs exists but no library analysis - run library analysis
    # If Ranger but not Bender outputs exists - run CellBender
    if(ranger.output.exists & bender.output.exists & !analysis.output.exists) {
      if(grepl("MAP", lib) | demux.output.exists) {
        analysis.run <- c(analysis.run, lib)  
      } else {
        missing.demux <- c(missing.demux, lib)
      }
    } else { if(ranger.output.exists & !bender.output.exists) {
      bender.run <- c(bender.run, lib)
    }}
  }
}
openxlsx::write.xlsx(details, args$job.details, rowNames=T)

if(length(setdiff(libs, rownames(details))) > 0) {
  message(paste("The following folders are not managed by", args$job.details))
  message(paste(setdiff(libs, rownames(details)), collapse = "\n"))
}

if(length(setdiff(rownames(details), libs)) > 0) {
  message(paste("No folder matching names of the following records in", args$job.details))
  message(paste(setdiff(rownames(details), libs), collapse = "\n"))
}

if(length(missing.demux) > 0) {
  message("The following folders are missing demultiplexing results:")
  message(paste(missing.demux, collapse = "\n"))
}

# -------------------------------------------------- #
# Create script files for jobs                       #
# -------------------------------------------------- #
if(length(bender.run) == 0) {
  message("No libraries need running CellBender")
} else {
  message(ifelse(args$execute, "\n\nExecuting CellBender on libraries:", "\n\nShould run CellBender on the following libraties:"))
  for (lib in bender.run) {
    if(args$execute) {
      if(file.exists(file.path(args$path, lib, "run.cellbender.sh"))) file.remove(file.path(args$path, lib, "run.cellbender.sh"))
      lapply(c("#!/bin/bash",
               "#SBATCH --partition=gpu.q",
               "#SBATCH --nodes=1",
               "#SBATCH --ntasks=1",
               "#SBATCH --gres=gpu:1",
               paste0("#SBATCH --mem=", details[lib,]$bender.ram, "GB"),
               paste0("#SBATCH --output=",file.path(args$path, lib,"%x.o%j")),
               paste0("rm -rf ", file.path(args$path, lib, args$bender.folder)),
               paste0("mkdir ", file.path(args$path, lib, args$bender.folder)),
               "source ~/miniconda/bin/activate CellBender",

               paste0('cellbender remove-background',
                      ' --input "',  file.path(args$path, lib, 'raw_feature_bc_matrix.h5'), '"',
                      ' --output "', file.path(args$path, lib, args$bender.folder, paste0(lib, '.h5')), '"',
                      ' --cuda --epochs 300 ',
                      ifelse(is.na(details[lib,]$expected.cells), '',
                             paste0('--expected-cells ', details[lib,]$expected.cells)),
                      paste0(' --total-droplets-included ', details[lib,]$total.droplets),
                      ' --learning-rate 0.00001 --z-dim 50')
      ),
      write, file.path(args$path, lib, "run.cellbender.sh"), append=T, ncolumns=1000)
      submission <- system(paste0("sbatch ", file.path(args$path, lib, "run.cellbender.sh")), intern = T)
      message(paste0(lib, ":\t", submission))
      details[lib,]$job.id <- gsub(".* ", "", submission)
      details[lib,]$job.name <- "CellBender"
      details[lib,]$job.status <- "Running"
      openxlsx::write.xlsx(details, args$job.details, rowNames=T)
    } else
      message(paste0("sbatch ", file.path(args$path, lib, "run.cellbender.sh")))
  }
}


if(length(analysis.run) == 0) {
  message("No libraries need running initial analysis")
} else {
  message(ifelse(args$execute, "\n\nExecuting library analysis on libraries:", "\n\nShould run library analysis on the following libraties:"))
  for (lib in analysis.run) {
    if(args$execute) {
      if(file.exists(file.path(args$path, lib, "run.library.analysis.sh"))) file.remove(file.path(args$path, lib, "run.library.analysis.sh"))
      lapply(c("#!/bin/bash",
               "#SBATCH --partition=elsc.q",
               "#SBATCH --nodes=1",
               "#SBATCH --ntasks=1",
               "#SBATCH --cpus-per-task=1",
               paste0("#SBATCH --mem=", details[lib,]$analysis.ram, "GB"),
               paste0("#SBATCH --output=",file.path(args$path, lib,"%x.o%j")),
               paste0("rm -rf ", file.path(args$path, lib, args$analysis.folder)),
               paste0("mkdir ", file.path(args$path, lib, args$analysis.folder)),
               paste0('Rscript ~/500/all/src/library.analysis.R',
                      ' -n "', lib, '"',
                      ' -l "', file.path(args$path, lib, args$bender.folder, paste0(lib, '_filtered.h5')),'"',
                      ' -d "', file.path(args$path, lib), '"',
                      ' -o "', file.path(args$path, lib, args$analysis.folder), '"',
                      ' -c "', details[lib,]$remove.lowqual.clusters, '"')),
             write, file.path(args$path, lib, "run.library.analysis.sh"), append=T, ncolumns=1000)
      submission <- system(paste0("sbatch ", file.path(args$path, lib, "run.library.analysis.sh")), intern = T)
      message(paste0(lib, ":\t", submission))
      details[lib,]$job.id <- gsub(".* ", "", submission)
      details[lib,]$job.name <- "LibraryAnalysis"
      details[lib,]$job.status <- "Running"
      openxlsx::write.xlsx(details, args$job.details, rowNames=T)
    } else
      message(paste0("sbatch ", file.path(args$path, lib, "run.library.analysis.sh")))
  }
}
