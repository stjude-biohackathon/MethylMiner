#!/bin/env Rscript
#########################################################################
# St. Jude BioHackathon 2022
# Team 2
#########################################################################

###############################################################################
# Load dependencies and helper scripts
.libPaths(c(.libPaths(),"./Rlib"))
print(.libPaths())
r = getOption("repos")
r["CRAN"] = "https://cloud.r-project.org/" #Use the cloud repo.
options(repos = r)

## Check if funr is installed (it is essential to get the source file path).
if(!'argparse' %in% .packages(all.available = TRUE)){
  message("The essential package 'argparse' is missing, installing it ...")
  install.packages('argparse')
}


if(!'funr' %in% .packages(all.available = TRUE)){
  message("The essential package 'funr' is missing, installing it ...")
  install.packages('funr')
}

## Get curent working directory and load scripts
suppressPackageStartupMessages(require(funr))
rootDir = dirname(funr::sys.script())
sourceDir = paste0(rootDir,"/R")
tmp =lapply( list.files(sourceDir,pattern = ".R",full.names = TRUE) , source)

## check if all the packages are installed
cli::cli_blockquote("Methylation QC and normalization")

cli::cli_h1("Loading and installing missing packages")
checkPackages()

###############################################################################
# Parse parameters
cli::cli_h1("Parse parameters and folder prep")
config  = parseParameters()
config$scriptsDir = rootDir

###############################################################################
# Prepare folders
res = prepareFolders(config)
config$outputDir =  res$outputDir
all_targets = res$all_targets



###############################################################################
# Load manifest
cli::cli_h1("Loading manifest")
manifest  = loadSourceFiles(config)

###############################################################################
# load data into minfi
cli::cli_h1("Create Minfi")
RGSet.all = LoadIntoMinfi(config, all_targets)

###############################################################################
# pre-Normalization QC
cli::cli_h1("Pre-normalization QC")
preNormalizationAndQC(config,RGSet.all)

#################################################################
# Quality control and normalization: sex check
cli::cli_h1("Sex check")
phenoData.all = sexCheck(config, RGSet.all)

#################################################################
# Quality control and normalization: Normalization
cli::cli_h1("Normalization")
target.MSet.norm.all = Normalization(config, RGSet.all)

#################################################################
# Quality control and normalization: Filter - separate autosomes from sex chromosomes and get beta values
cli::cli_h1("Separating automsomes")
seprateAutosomes(config,RGSet.all, target.MSet.norm.all, manifest, phenoData.all)

logger::log_info("Step1: QC & Normalization [DONE!]")




