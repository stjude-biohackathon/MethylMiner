#!/usr/bin/r Rscript
# Standard Pipeline for analysis of methylation array data
# Part 1: Quality Control and Normalization using minfiR
# Christy LaFlamme
# 02-28-22
###############################################################################
# command: Rscript 1_methylation_QC_norm.R -w workDir -n name_of_run
# arguments:
# workDir - working directory should be main folder containing the following folders:
# data directory should contain all raw .idats for analysis and metadata csv file
# output directory to direct output (if folder does not exist, the script will create one)
# within the output directory, a folder will be created for this QC/norm output detailing the number of samples
# name_of_run - name of QC/norm run

# args = commandArgs(trailingOnly=TRUE) # get arguments

# default normalization method is preprocessIllumina
# if other normalization is preferred, must change manually
###############################################################################
library(optparse)

option_list = list(
  make_option(c("-w", "--workDir"), type="character", default=NULL, 
              help="working directory file path", metavar="character"),
  make_option(c("-n", "--runName"), type="character", default="run", 
              help="output folder name", metavar="character"),
  make_option(c("-m", "--manifest"), type="character", default="run", 
              help="manifest file", metavar="character")

); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$workDir)){
  print_help(opt_parser)
  stop("Working directory must be set.", call.=FALSE)
}

###############################################################################
# Load dependencies
library(dplyr)
library(tidyverse)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(data.table)
library(readxl)

# Load manifest

manifest = toString(opt$manifest)
EPIC.hg19.manifest = read.table(manifest, sep="\t") 
colnames(EPIC.hg19.manifest) <- EPIC.hg19.manifest[1,] # rename the column names as the first row
EPIC.hg19.manifest <- EPIC.hg19.manifest[-1,] # remove the first row
rownames(EPIC.hg19.manifest) <- seq(length=nrow(EPIC.hg19.manifest)) # restart row numbering

#source("/Users/claflamm/Downloads/source_files.R")

  #################################################################
  # set directories
  workDir <- toString(opt$workDir)

  dataDir <- paste(workDir, "/", sep = "") 
  print(dataDir)
  
  all_targets <- read.metharray.sheet(dataDir) # read sample sheet containing metadata
  
  outputDir <- paste(getwd(), "/", sep = "") # if output directory does not exist, create one
  
  if (!file.exists(outputDir)) {
    dir.create(outputDir)
  }
  
  name <- toString(opt$runName)
  
  QC_run_folder <- paste(outputDir, name, "_preprocessIllumina/", sep = "") # make specific folder for this QC run
  
  if (!file.exists(QC_run_folder)) {
    dir.create(QC_run_folder)
  }
  
  outputDir <- QC_run_folder # set new output directory as this folder

  #################################################################
  # Quality control and normalization: load data into minfi
  
  RGSet.all <- read.metharray.exp(targets = all_targets, force=TRUE) # get RGChannelSet object
  
  if (!file.exists(paste(outputDir, "raw_data", sep = ""))) { # folder for output using raw data
    dir.create(paste(outputDir, "raw_data", sep = ""))
  }
  
  saveRDS(RGSet.all, file = paste(outputDir,"raw_data/RGSet.all.RDS", sep = "")) # save data to RDS
  
  annotation(RGSet.all) # check annotation
  
  phenoData.all <- pData(RGSet.all) # get the phenotype data
  phenoData.all <- as.data.frame(phenoData.all)

  saveRDS(phenoData.all, file = paste(outputDir,"phenoData.all.RDS", sep = "")) # save data to RDS
  
  #################################################################
  # Quality control and normalization: pre-Normalization QC
  
  MSet.raw.all <- preprocessRaw(RGSet.all) # get raw methyl set
  
  qc.raw.all <- getQC(MSet.raw.all) # get quality control - methylated median and unmethylated median
  
  # get quality control plot
  pdf(paste(outputDir,"raw_data/plotQC_raw.pdf", sep = ""))
  plotQC(qc.raw.all)
  dev.off()
  
  # NOTE: should have step in here that removes failed samples; however, currently not sure how to tell which samples failed except based off of the picture with the indexes
  
  # plot density plots
  pdf(paste(outputDir,"raw_data/density_beta_raw.pdf", sep = ""))
  densityPlot(MSet.raw.all, sampGroups = phenoData.all$Sample_Group)
  dev.off()

  # generate quality control report
  qcReport(RGSet.all, pdf= paste(outputDir,"raw_data/qcReport.pdf", sep = ""))
  
  #################################################################
  # Quality control and normalization: sex check
  
  GRSet <- mapToGenome(RGSet.all) # get genomic ratio set 
   
  if (!file.exists(paste(outputDir, "sex_prediction", sep = ""))) { # folder for output using raw data
    dir.create(paste(outputDir, "sex_prediction", sep = ""))
  }
  
  predictedSex <- getSex(GRSet, cutoff = -2) # get the predicted sexes
  
  GRSet <- addSex(GRSet, sex = predictedSex) # add the sex prediction info to the genomic ratio set
  
  pdf(paste(outputDir,"sex_prediction/sexplot.pdf", sep = "")) # plot the predicted sex
  plotSex(GRSet)
  dev.off()
  
  # add predicted sex information to phenoData and check to see if reported sex is equal to predicted sex
  phenoData.all$Predicted_Sex <- as.factor(predictedSex$predictedSex)
  phenoData.all$Predicted_Sex <- factor(phenoData.all$Predicted_Sex, levels = c("F", "M"),labels = c("Female","Male")) # change the factor names
  phenoData.all$Sex_Concordance <- ifelse(as.character(phenoData.all$Predicted_Sex) == as.character(phenoData.all$Reported_Sex), TRUE, FALSE)
  
  # save sex predictions to text file
  write.table(phenoData.all[,c("Sample_Name", "Sample_ID", "Reported_Sex", "Predicted_Sex", "Sex_Concordance")], file = paste(outputDir, "sex_prediction/sexPredictionsComparison.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  
  # save sex mismatch samples to text file
  discordant.sex <- phenoData.all %>% filter(Sex_Concordance == FALSE)
  write.table(as.data.frame(discordant.sex), file = paste(outputDir, "sex_prediction/sexPredictionsComparison.discordant.txt", sep = ""),  sep = "\t", quote = FALSE, row.names = FALSE) 
  
  #################################################################
  # Quality control and normalization: Normalization
  
  # background subtraction and control normalization
  MSet.norm.all <- preprocessIllumina(RGSet.all, bg.correct = TRUE, normalize = "controls") 
  
  n_samples <- ncol(MSet.norm.all)
  n_probes <- nrow(MSet.norm.all)
  
  # A detection p-value is returned for every genomic position in every sample. 
  detect.levels.all <- detectionP(RGSet.all, type="m+u") # The m+u method compares the total DNA signal (Methylated + Unmethylated) for each position to the background signal level. 
  
  # designate failed samples 
  failed <- detect.levels.all > 0.01 
  fraction_failed <- colMeans(failed) # fraction of failed positions per sample (# true/# total)
  fraction_failed_50 <- sum(rowMeans(failed) > 0.5) # how many positions failed in >50% of samples?
  print(paste(fraction_failed_50, "out of", n_probes, "probes failed in >50% of the samples"))
  
  if (!file.exists(paste(outputDir, "normalized_data", sep = ""))) { # folder for output using raw data
    dir.create(paste(outputDir, "normalized_data", sep = ""))
  }
  
  # plot fraction of failed positions per sample
  pdf(paste(outputDir,"normalized_data/Fraction_Failed_Positions_PerSample.pdf", sep = ""))
  barplot(sort(fraction_failed), main="EPIC Fraction of Failed Positions per Sample", ylab="Fraction Failed Positions")
  dev.off()
  
  # how many samples have significant p values at each site?
  detect.sd.all <- rowSums(detect.levels.all < 0.01) # sums the number of samples for which that probe[row] has a "significant" p value
  # the number associated with each probe is equivalent to the number of samples with significant p values at that site
  print(paste("The average number of samples with significant p-values across all probes is", mean(detect.sd.all), "out of", n_samples))
  
  # which sites have >90% success rate of samples?
  detect.ind.sd.all <- which(detect.sd.all > (0.9*n_samples)) # determine which rows have > 90% of the total number of samples with "significant" p-values
  
  print(paste(length(detect.ind.sd.all), "out of", n_probes, "probes have >90% of the total number of samples with significant p-values"))
  
  detect.ind.sd.all <- as.matrix(detect.ind.sd.all) # convert to matrix
  
  # determine the overlap between the normalized methyl set and the desired sites
  overlap.all <- intersect(rownames(MSet.norm.all), rownames(detect.ind.sd.all))
  # length(overlap.all) == length(detect.ind.sd.all) # should be same length as detect.ind.sd.all
  
  target.MSet.norm.all <- MSet.norm.all[overlap.all,] # create a new methyl set containing only the desired sites (>90% success rate)
  
  number.sites.dropped <- nrow(MSet.norm.all) - nrow(target.MSet.norm.all) # check the number of probes dropped
  
  # print(paste("Dropping", number.sites.dropped, "probes..."))
  print(paste("Of", nrow(MSet.norm.all), "total probes, taking", nrow(target.MSet.norm.all), "probes with significant p-values in >90% of samples"))
  
  target.ratioSet.all <- ratioConvert(target.MSet.norm.all, what = "both", keepCN = TRUE) # convert to a ratio set

  GSet.all <- mapToGenome(target.ratioSet.all) # map the target ratio set to the genome
  GSet.all <- dropLociWithSnps(GSet.all) # remove loci with snps
  
  saveRDS(GSet.all, paste(outputDir,"normalized_data/GSet.all.RDS", sep = "")) # save data to RDS 
  
  #################################################################
  # Quality control and normalization: Filter - separate autosomes from sex chromosomes and get beta values
  
  # since the cohort has both males and females, separate the autosomes and perform sexchr analysis separately on males and femaless
  autosomes <- !(featureNames(target.MSet.norm.all) %in% EPIC.hg19.manifest$probeID[EPIC.hg19.manifest$CpG_chrm %in% c("chrX", "chrY")]) # ! = not; keep everything but sex chromosomes
  sexchr <- (featureNames(target.MSet.norm.all) %in% EPIC.hg19.manifest$probeID[EPIC.hg19.manifest$CpG_chrm %in% c("chrX", "chrY")])
  
  MSet.norm.clean.auto <- target.MSet.norm.all[autosomes,] # autosomes for males and females
  
  phenoData.female <- phenoData.all[phenoData.all$Predicted_Sex == "Female",]
  females <- row.names(phenoData.female) 
  MSet.norm.clean.females.sexchr <- target.MSet.norm.all[sexchr,females] 
  
  phenoData.male <- phenoData.all[phenoData.all$Predicted_Sex == "Male",]
  males <- row.names(phenoData.male)
  MSet.norm.clean.males.sexchr <- target.MSet.norm.all[sexchr,males] 
  
  saveBeta <- function(MSet, phenoData, filename) {
    
    #################################################################################################################
    # get beta values
    GSet <- mapToGenome(MSet) # convert to genomic methyl set
    GSet <- dropLociWithSnps(GSet) # remove loci with snps
    beta <- getBeta(GSet) # get beta values
    # saveRDS(beta, file = paste(outputDir, "normalized_data/beta.", filename, ".RDS", sep = "")) # save data to RDS
    
    #################################################################################################################
    # plot a density plot of the beta values
    pdf(paste(outputDir, "normalized_data/beta.density.bysample.", filename, ".pdf", sep = ""))
    densityPlot(beta, sampGroups=phenoData$Sample_Name)
    dev.off()
    
    pdf(paste(outputDir, "normalized_data/beta.density.bygroup.", filename, ".pdf", sep = ""))
    densityPlot(beta, sampGroups=phenoData$Sample_Group)
    dev.off()
    
    #################################################################################################################
    # annotate beta values
    beta <- as.data.frame(cbind(rownames(beta), beta)) # add probeID info as first column
    colnames(beta)[1] <- "probeID"
    
    filtered.probes <- beta$probeID # filter the EPIC manifest data to only contain the filtered probes
    
    # filter the manifest data to only contain the filtered probes
    filtered.probes.manifest <- EPIC.hg19.manifest[EPIC.hg19.manifest$probeID %in% filtered.probes,] # designate the EPIC manifest data for the filtered probes
    filtered.probes.manifest <- filtered.probes.manifest %>% select(5, 1:3)
    
    # merge the objects at the probe IDs
    beta <- merge(filtered.probes.manifest, beta, by.x = "probeID", by.y = "probeID", all.x = FALSE, all.y = FALSE)
    write.table(beta, file = paste(outputDir, "normalized_data/", filename, ".beta.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t") # save beta values as txt file
    
  }
  
  saveBeta(MSet.norm.clean.auto, phenoData.all, "autosomes")
  saveBeta(MSet.norm.clean.females.sexchr, phenoData.female, "female.sexchr")
  saveBeta(MSet.norm.clean.males.sexchr, phenoData.male, "male.sexchr")
  






























