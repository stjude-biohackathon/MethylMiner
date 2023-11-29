# Added features to DMR analysis for MethFlow Pipeline
# 04-25-22
# updates 08-10-22
# turned into functions on 08-15-22
# Christy LaFlamme

#########################################################
# load dependencies
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggprism)
library(GenomicRanges)
library(reshape2)

#########################COORDINATES ANNOTATION#################################
add_coordiantes <- function(dmrlist) {
  
  # add useful information about the coordinates
  dmrlist$DMR_coordinates <- paste(dmrlist$seqnames, ":", dmrlist$start, "-", dmrlist$end, sep = "") # add column with location for convenient copy paste
  dmrlist$DMR_coordinates_dir <- paste(dmrlist$seqnames, ":", dmrlist$start, "-", dmrlist$end, "[", dmrlist$direction, "]",sep = "") # add column for location+direction defining "unique" DMRs
  dmrlist$DMR_order <- 1:nrow(dmrlist)
  
  return(dmrlist)
}

#########################CTCF ANNOTATION########################################
annotate_CTCF <- function(dmrlist, annotDir) {
  
  CTCF.overlap <- as.data.frame(fread(paste(annotDir, "EPIC.hg19.CTCF.overlap.tsv", sep = ""),)) # load the CTCF annotation
  
  dmrlist$CTCF_overlap <- NA 
  # add information about overlapping CTCF region
  for (i in 1:nrow(dmrlist)) { # loop over the dmrlist
    for (j in 1:nrow(CTCF.overlap)) { # loop over the rows of the CTCF file
      if(dmrlist$seqnames[i] == CTCF.overlap$CpG_chrm[j]) { # if the seqnamesosomes are equal
        if (dmrlist$start[i] <= CTCF.overlap$CpG_beg[j] & dmrlist$end[i] >= CTCF.overlap$CpG_end[j]) { # make logical comparison: if the start of the dmr is less than or equal to the start of the CTCF boundary AND if the end of the DMR is greater than or equal to the end of the CTCF boundary, then add the CTCF boundary to the annotation information
          dmrlist$CTCF_overlap[i] <- CTCF.overlap$CTCF[j]
        }
      }
    }
  }
  
  # add information about number of overlapping CTCF probes
  colnames(CTCF.overlap)[1:3] = c("seqnames","start","end")
  CTCF.overlap <- GRanges(CTCF.overlap) # convert to GRanges object
  dmrlist <- GRanges(dmrlist) # convert to GRanges object
  
  ovp <- findOverlaps(dmrlist, CTCF.overlap) # find overlapping regions
  dmr.hits = table(queryHits(ovp)) # get a table of the number of overlapping regions for each dmr
  
  dmrlist$CTCF_overlap_nprobes = 0 # add column to store number of overlapping probes
  dmrlist$CTCF_overlap_nprobes[as.numeric(names(dmr.hits))] = dmr.hits # add this information from the table
  
  return(as.data.frame(dmrlist))
}

#########################450K ANNOTATION########################################
# 450K CONTROL DATA ANNOTATION
# add n probes and 450K probe info (whether probes were included on the 450k array

annotate_450k <- function(dmrlist, data, annotDir) {
  anno450kmanifest <- fread(paste(annotDir, "HM450.hg19.manifest.tsv", sep = ""))
  
  dmrlist <- as.data.frame(dmrlist)
  dmrlist$probeIDs <- NA
  dmrlist$nprobes <- NA
  dmrlist$check.probes.450k <- NA 
  
  for (i in 1:nrow(dmrlist)){ # write a for loop to check all if probes are included 450k array
    chr = NULL
    start = NULL
    end = NULL
    probes= NULL
    n.probes = NULL
    probe.logical = NULL
    
    chr <- dmrlist$seqnames[i] # set the chromosome
    
    start <- which(data[data$CpG_chrm == chr,]$CpG_beg == paste(dmrlist$start[i])) # get the epivariation data index for the dmr
    end <- which(data[data$CpG_chrm == chr,]$CpG_end == paste(dmrlist$end[i]))
    
    probes <- data[start:end,probeID] # get the probe IDs for the dmr
    dmrlist$probeIDs[i] <- toString(probes)
    
    n.probes <- length(probes) # count the number of probes
    dmrlist$nprobes[i] <- n.probes
    
    probe.logical <- probes %in% anno450kmanifest$probeID # get a logical statement indicating whether or not probe is included on 450k array
    dmrlist$check.probes.450k[i] <- toString(probe.logical) # append the logical to the dataframe as a string
    
  }
  return(dmrlist)
}

#########################BRAIN EXP ANNOTATION########################################

annotate_brain_exp <- function(dmrlist, annotDir) {
  
  gene.brain.exp <- fread(paste(annotDir, "rna_brain_gtex.tsv", sep = "")) # read in RNA brain GTex data
  gene.brain.exp <- gene.brain.exp[,c("Gene name", "nTPM")] # get the gene name and normalized transcripts per million column
  colnames(gene.brain.exp)[2] <- "Brain_Gene_Expression_nTPM"
  gene.brain.exp <- as.data.frame(gene.brain.exp[, lapply(.SD, mean), by= "Gene name"]) # take the average across brain regions
  
  dmrlist <- merge(dmrlist, gene.brain.exp, by.x = "Gene.Name", by.y = "Gene name", all.x = TRUE, all.y = FALSE, sort = FALSE) 
  
}


#########################OMIM ANNOTATION########################################
annotate_OMIM <- function(dmrlist, annotDir) {
  
  OMIM <- fread(paste(annotDir, "genemap2.txt", sep = ""), skip = 3)
  OMIM <- as.data.frame(OMIM %>% separate_rows(`Gene Symbols`))
  
  dmrlist <- merge(dmrlist, OMIM[c("Gene Symbols", "MIM Number", "Phenotypes")], by.x = "Gene.Name", by.y = "Gene Symbols", all.x = TRUE, all.y = FALSE, sort = FALSE) # note that 'Gene Name' was changed to 'Gene.Name' by GRanges in annotation_CTCF function
  
  return(dmrlist)
}

#########################ADD PHENO########################################
add_pheno <- function(dmrlist, phenoData.all) {
  dmrlist.pheno <- merge(dmrlist, phenoData.all, by.x = "EpiInd", by.y = "Sample_ID", all.x = TRUE, all.y = FALSE, sort = FALSE)
  dmrlist.pheno <- dmrlist.pheno[order(dmrlist.pheno$DMR_order),]
  
  return(dmrlist.pheno)
}

#########################SAMPLE SUMMARY########################################
sample_summary <- function(dmrlist, phenoData.all, workDir, dmrDir, chromosome) {
  
  ########################################################
  # get summary of DMRs by the sample
  sampleID.summary <- as.data.frame(sort(table(dmrlist$EpiInd), decreasing = TRUE))
  sampleID.summary <- merge(sampleID.summary, phenoData.all[,c("Sample_ID", "Sample_Name", "Sample_Type", "Sample_Plate")], by.x = "Var1", by.y = "Sample_ID", all.x = TRUE, all.y = FALSE, sort = FALSE)
  colnames(sampleID.summary)[which(names(sampleID.summary) == "Var1")] <- "Sample_ID"
  colnames(sampleID.summary)[which(names(sampleID.summary) == "Freq")] <- "nDMRsperSample"
  write.table(sampleID.summary, file = paste(workDir, dmrDir, chromosome, "_sample_summary", "/sample.summary.txt", sep = ""), quote = F, sep = "\t", row.names = F)
  
  return(sampleID.summary)
}

#########################PLOTTING DISTRIBUTION########################################
plot_distribution <- function(workDir, dmrDir, chromosome) {
  
  sampleID.summary <- fread(paste(workDir, dmrDir, chromosome, "_sample_summary", "/sample.summary.txt", sep = ""))
  
  # plot distribution of DMRs across samples
  pdf(paste(workDir, dmrDir, chromosome, "_sample_summary", "/DMR.distribution.pdf", sep = ""))
  
  # plot distribution of DMRs across samples
  g <- ggplot(data=as.data.frame(sampleID.summary), aes(x = nDMRsperSample ))+ # plot the data # Freq = number of DMRs
    geom_histogram(binwidth = 10)+
    expand_limits(y = 450)+
    ggtitle("DMR Distribution")+
    xlab("Number of DMRs")+
    ylab("Number of Patients")+ # y axis label
    theme_prism(base_size = 20)+ # graphpad prism themed
    scale_y_continuous(expand = c(0, 0))+ # make the zero on the y axis match the x axis
    scale_x_continuous(breaks = c(1, 100, 1000, 5000))+
    stat_bin(aes(y=..count.., label=..count..), geom="text", vdmrust=-0.4, size = 6, binwidth = 30)
  
  print(g)
  
  dev.off()
}

#########################ANNOTATE NUM DMRs########################################
annotate_nDMRs <- function(dmrlist.pheno, sampleID.summary, workDir, dmrDir, dmrFile) {
  # add number of DMRs per patient to annotated dmrlist
  dmrlist.pheno <- merge(dmrlist.pheno, sampleID.summary[,c("Sample_ID", "nDMRsperSample")], by.x = "EpiInd", by.y = "Sample_ID", all.x = TRUE, all.y = TRUE)
  dmrlist.pheno <- dmrlist.pheno[order(dmrlist.pheno$DMR_order),]
  dmrFilePath <- paste(workDir, dmrDir, "annotations/", dmrFile, sep = "")
  write.table(dmrlist.pheno, file = sub(".tsv", ".meta.tsv", dmrFilePath), quote = F, sep = "\t", row.names = F)
  
  return(dmrlist.pheno)
}

#########################SAMP100 ########################################
samp100 <- function(sampleID.summary, workDir, dmrDir, chromosome) {
  # select samples that display over 100 DMRs
  samp100 <- as.data.frame(sampleID.summary %>% filter(nDMRsperSample > 100) %>% select(c(Sample_ID, Sample_Name, nDMRsperSample)))
  sentrix100 <- as.data.frame(as.character(samp100$Sample_ID))
  colnames(sentrix100)[1] <- "Sample_ID"
  write.table(sentrix100, paste(workDir, dmrDir, chromosome, "_sample_summary", "/samples.with.over.100.DMRs.txt", sep = ""), quote = F, sep = "\t", row.names = F, col.names = F)
}

#########################CONSOLIDATE DMRs########################################
consolidate_dmrs <- function(dmrlist.pheno, workDir, dmrDir, chromosome) {
  
  # consolidate DMR list (by chr, start, stop, and direction) for DMR plots so that each region of differential methylation is plotted only once
  dmrlist.consolidate <- as.data.table(dmrlist.pheno)[, toString(Sample_Name), by = list(seqnames, start, end, direction, DMR_coordinates_dir, DMR_size, Annotation, `Distance.to.TSS`, `Gene.Name`)]
  colnames(dmrlist.consolidate)[which(names(dmrlist.consolidate) == "V1")] <- "Sample_Name"
  dmrlist.consolidate$DMR_number <- 1:nrow(dmrlist.consolidate)
  write.table(dmrlist.consolidate, file = paste(workDir, dmrDir, chromosome, "_dmr_plots",  "/dmrlist.plot.reference.txt", sep = ""), quote = F, sep = "\t", row.names = F)
  # dmrlist.consolidate <- fread(paste(workDir, dmrDir, "autosomes/dmr_plots/dmrlist.plot.reference.txt", sep = ""))
  
  # add DMR number information to annotated dmrlist
  dmrlist.pheno <- merge(dmrlist.pheno, dmrlist.consolidate[,c("DMR_coordinates_dir","DMR_number")], by.x = "DMR_coordinates_dir", by.y = "DMR_coordinates_dir", all.x = TRUE, all.y = FALSE, sort = FALSE)
  dmrlist.pheno <- dmrlist.pheno[order(dmrlist.pheno$DMR_order),]
  dmrFilePath <- paste(workDir, dmrDir, "annotations/", dmrFile, sep = "")
  write.table(dmrlist.pheno, file = sub(".tsv", ".plotref.tsv", dmrFilePath), quote = F, sep = "\t", row.names = F)
  
  return(dmrlist.consolidate)
}

#########################PLOTTING SCRIPT########################################
# DMR plotting script (Paras - Sharp lab)

plotDMR <- function(data, samples, startPos, endPos, chrCol, startCol, title = ""){
  data = data %>% slice(startPos:endPos)
  epivarSamples = data %>% select(Sign_individuals_t0.1_n3_w1k_BOTH) %>% 
    filter(!is.na(Sign_individuals_t0.1_n3_w1k_BOTH)) %>% # filter out N/A
    separate_rows(Sign_individuals_t0.1_n3_w1k_BOTH, sep = ",") %>% # generate separate rows for comma separated values
    distinct %>% pull
  
  data = data %>% 
    melt(id.vars = setdiff(colnames(data), samples), variable.name = "sampleID", value.name = "Beta") %>%
    mutate(status = ifelse(sampleID %in% epivarSamples, "DMR-carrier", "normal"))
  
  g = ggplot(data =  data, aes_string(x=startCol, y="Beta")) +
    geom_line(aes(group=sampleID,size=status, linetype = status, colour = status)) +
    scale_colour_manual(values = c("DMR-carrier" = "#cf3a36","normal" = "black")) + # #cf3a36 used to be "red"
    scale_size_manual(values = c("DMR-carrier"=0.9, "normal" = 0.2 )) +
    scale_linetype_manual(values = c("DMR-carrier"="solid", "normal" = "dashed" )) +
    coord_cartesian(ylim = c(0,1)) +
    xlab(paste0("Genomic Position at ",data[1,chrCol])) +
    ylab("Methylation Value") +
    labs(title = title) +
    theme_bw()
  g
}

#########################PLOTTING LOOP########################################
# loop for plotting DMR plots

loopDMR <- function(dmrlist, data, samples, workDir, dmrDir, chromosome) {
  # plotting loop
  for (dmr in 1:nrow(dmrlist)) {
    
    num = NULL
    dir = NULL
    chrom = NULL
    start = NULL
    end = NULL
    width = NULL
    gene = NULL
    annot = NULL
    samp = NULL
    
    num = dmrlist$DMR_number[dmr]
    dir = dmrlist$direction[dmr]
    chrom = dmrlist$seqnames[dmr]
    start = dmrlist$start[dmr]
    end = dmrlist$end[dmr]
    width = dmrlist$DMR_size[dmr]
    gene = gsub("/", "_", dmrlist$`Gene.Name`[dmr]) # replace any "/" with _ or the path will split the name of the file into a directory
    annot = dmrlist$`Annotation`[dmr]
    samp = dmrlist$`Sample_Name`[dmr]
    
    file = paste(workDir, dmrDir, chromosome, "_dmr_plots", "/by_dmr/", num, "_", dir, "_", chrom, "_", start, "_", end, "_", width, "_", gene, "_", annot, "_", samp, ".pdf", sep = "")
    
    if(!file.exists(file)) {
      
      Start_DMR = NULL
      End_DMR = NULL
      
      Start_DMR = which(data$CpG_beg == paste(dmrlist[dmr, "start"])) 
      End_DMR = which(data$CpG_end == paste(dmrlist[dmr, "end"]))
      
      pdf(file)
      
      DMRplot <- plotDMR(data, samples, Start_DMR, End_DMR, "CpG_chrm","CpG_beg", title = gene)
      print(DMRplot)
      
      dev.off()
      print(file)
      
    }
    
  }
  
}

#########################################################
# write a function for dmr analysis

dmr_analysis <- function(workDir, dmrDir, annotDir, dmrFile, epiFile, chromosome, n, n_subset) {
  
  #########################################################
  # load dependencies
  library(data.table)
  library(dplyr)
  library(tidyr)
  # library(tidyverse)
  library(ggplot2)
  # library(ggprism)
  library(GenomicRanges)
  library(reshape2)
  
  #########################################################
  # load data
  
  print("Loading phenoData and dmrlist")
  # phenoData.all <- readRDS(paste(workDir, "phenoData.all.RDS", sep = "")) # meta data
  phenoData.all <- readRDS("/home/claflamm/methylation/raw_idats/1311_all_metharray_FINAL_cohort_run/output/1311_all_metharray_FINAL_cohort_run/1311_all_metharray_FINAL_cohort_run_preprocessIllumina/normalized_data/betaValues/phenoData.all.RDS") # adjustment for pl22-24 analysis
  dmrlist <- fread(paste(workDir, dmrDir, "annotations/", dmrFile, sep = "")) # annotated dmr list
  print("Loaded")
  #########################################################
  # add useful information about coordiantes
  dmrlist<- add_coordiantes(dmrlist)
  print("DMR coordinates added")
  #########################################################
  # correct the CTCF annotation
  print("Starting CTCF annotation")
  dmrlist <- annotate_CTCF(dmrlist, annotDir)
  print("CTCF annotation complete")
  #########################################################
  # load epivariation data (beta value matrix)
  # this is required for 450k annotation
  print("Loading epivariation data")
  data <- fread(paste(workDir, dmrDir, epiFile, sep = ""))
  print("Loaded")
  #########################################################
  # add 450k annotation
  print("Starting 450k annotation")
  dmrlist <- annotate_450k(dmrlist, data, annotDir)
  print("450k annotation complete")
  #########################################################
  # add gene expression information
  # Ensembl gene identifier ("Gene"), analysed sample ("Tissue"), transcripts per million ("TPM"), protein-transcripts per million ("pTPM") and normalized expression ("nTPM").
  print("Starting brain gene expression annotation")
  dmrlist <- annotate_brain_exp(dmrlist, annotDir)
  print("Brain gene expression annotation complete")
  #########################################################
  # add OMIM information for genes
  print("Starting OMIM annotation")
  dmrlist <- annotate_OMIM(dmrlist, annotDir)
  print("OMIM annotation complete")
  ########################################################
  # annotate the dmrlist with phenotype metadata (e.g. sample name) and number of DMRs per sample
  print("Starting phenoData annotation")
  dmrlist.pheno <- add_pheno(dmrlist, phenoData.all)
  print("phenoData annotation complete")
  #########################################################
  # generate sample summaries
  # create folder for sample summary information
  sample_summary_folder <- paste(workDir, dmrDir, chromosome, "_sample_summary",  sep = "")
  
  if (!file.exists(sample_summary_folder)) {
    dir.create(sample_summary_folder)
  }
  
  # static DMR distribution plot
  print("Starting sample summary")
  sampleID.summary <- sample_summary(dmrlist, phenoData.all, workDir, dmrDir, chromosome)
  plot_distribution(workDir, dmrDir, chromosome)
  
  # annotate dmrlist.pheno with nDMRs per sample
  dmrlist.pheno <- annotate_nDMRs(dmrlist.pheno, sampleID.summary, workDir, dmrDir, dmrFile)
  
  # get info for samples with over 100 DMRs
  samp100(sampleID.summary, workDir, dmrDir, chromosome)
  print("Sample summary complete")
  ########################################################
  # generate static DMR plots
  # create folder for dmr plot information
  dmr_plots <- paste(workDir, dmrDir, chromosome, "_dmr_plots", sep = "")
  
  if (!file.exists(dmr_plots)) {
    dir.create(dmr_plots)
  }
  
  # consolidate dmrs
  dmrlist <- consolidate_dmrs(dmrlist.pheno, workDir, dmrDir, chromosome)
  
  # get sample info
  samples <- colnames(data)
  samples <- tail(samples, -4) # remove probe, chr, start, stop
  samples <- head(samples, -5) # remove quantile and significance information
  
  # create folder for plotting the dmrs "by_dmr" as opposed to "by_sample" etc.
  by_dmr <- paste(workDir, dmrDir, chromosome, "_dmr_plots", "/by_dmr", sep = "")
  
  if (!file.exists(by_dmr)) {
    dir.create(by_dmr)
  }
  
  #########################################################
  # ERROR handling
  # detach("package:minfi")
  # detach("package:bumphunter")
  # detach("package:Biostrings")
  # detach("package:XVector")
  # detach("package:SummarizedExperiment")
  detach("package:GenomicRanges")
  detach("package:GenomeInfoDb")
  detach("package:IRanges")
  detach("package:S4Vectors")
  #########################################################
  # plot DMR plots in DMR loop
  print("Starting DMR plotting")
  loopDMR(dmrlist, data, samples, workDir, dmrDir, chromosome)
  print("DMR plotting complete")
}

##################USE THE FUNCTION#######################
# set variables to NULL
n = NULL
n_subset = NULL
workDir = NULL
dmrDir=NULL
annotDir= NULL
chromosome = NULL
window = NULL
stringency = NULL
dmrFile = NULL
epiFile = NULL

#########################################################
# set parameters based on analysis run
n = "1156" # file naming convention (optional)
n_subset = "1156_filtered_blood_PBMC_PCAoutliers_rm_CellTypeCompoutliers_rm_samp100DMR_rm" # file naming convention (optional)
window = 1000 # set during pipeline
max = "0.99" # set during pipeline
min = "0.01" # set during pipeline
chromosome = "autosomes" # autosomes or sexchr

# set working directories
# workDir = paste0("/home/claflamm/methylation/raw_idats/", n, "_all_metharray/output/", n, "_all_metharray_dev/", n, "_all_metharray_dev_preprocessIllumina/")
workDir = "/home/claflamm/methylation/raw_idats/1311_all_metharray_FINAL_cohort_run/output/1311_all_metharray_FINAL_cohort_run/1311_all_metharray_FINAL_cohort_run_preprocessIllumina/"
dmrDir = paste0("normalized_data/dmr_", window, "_", min, "_", max, "/")
annotDir <- "/home/claflamm/methylation/annotations/" # directory containing annotations (can be downloaded from github page)
################AUTOSOMES#####################
# set names of files
dmrFile = paste0(gsub("_", ".", chromosome), ".beta.txt.sorted_", window, "_", max, "_", min, "_findEpivariation.sig_dmr.anno.tsv")
epiFile = paste0(gsub("_", ".", chromosome), ".beta.txt.sorted_", window, "_", max, "_", min, "_findEpivariation.txt.sig")

dmr_analysis(workDir, dmrDir, annotDir, dmrFile, epiFile, chromosome, n, n_subset) # run command

