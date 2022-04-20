installMissingPkg <- function(...){

  requiredPackages<-unlist(list(...))

  if(!all(requiredPackages %in% .packages(all.available = TRUE))){
    missingPackages = setdiff(requiredPackages, .packages(all.available = TRUE))
    pak::pkg_install(missingPackages)
    missingPackages = setdiff(requiredPackages, .packages(all.available = TRUE))
    BiocManager::install(missingPackages,update=FALSE)
    missingPackages = setdiff(requiredPackages, .packages(all.available = TRUE))
    for(pkg in missingPackages){
      if( !is.element(pkg, .packages(all.available = TRUE)) ) {
        message(paste0(pkg," missing, installing it ..."))
        BiocManager::install(pkg,update=FALSE)
      }
    }    
  }

}


using<-function(..., character.only = FALSE) {

    if(character.only){
      libs <- unlist(list(...))
    }else{
      libs <- unlist(list(as.character(substitute(...) ) ))
    }
  
    req<-unlist( lapply(libs,function(p) suppressPackageStartupMessages(require(p,character.only=TRUE)) ) )
    need<-libs[req==FALSE]
    if(length(need)>0){ 
        installMissingPkg(need)
        lapply(need,require,character.only=TRUE)
    }
}


checkPackages <- function(){
  if(!suppressPackageStartupMessages(require("pak",character.only = TRUE))){
    install.packages("pak", repos = "https://r-lib.github.io/p/pak/dev/")
  }

  ## Add any required packages to this list
  requiredPackages <- c("cli","optparse","yaml","dplyr","tidyverse","glue","fs",
                        "minfi","IlluminaHumanMethylationEPICmanifest","processx",
                        "IlluminaHumanMethylationEPICanno.ilm10b4.hg19","openxlsx",
                        "GenomicRanges",
                        "data.table","readxl","ggplot2","logger")

  # just in case a package was written twice                    
  requiredPackages = unique(requiredPackages)   

  suppressPackageStartupMessages(using(requiredPackages, character.only=TRUE))  

  log_layout(layout_glue_colors)
  log_threshold(TRACE)
}

#' Function to read the configuration of the run
#'
#' The configuration file is a yaml fomated file, in which the parameters of each stage are
#' specified.
#' @return
#' A named list, each containing the parameters of the named stage
#' @export
#'
#' @examples
parseParameters <- function(){
    # Parse paramters
  parser <- argparse::ArgumentParser(
    description = "Main script to generated the initial QC and normalization."
  ); 

  parser$add_argument(
    "-w", "--workDir",
    dest='workDir',
    help = "Path to the working directory",
    type = "character"
  );

  parser$add_argument(
    "-n", "--runName",
    dest='runName',
    help = "Output folder name",
    type = "character"
  );

  parser$add_argument(
    "-g","--genome",
    dest='genome',
    help = "Genome build. Default: hg19",
    type = "character",
    default='hg19'
  );

  parser$add_argument(
    "-m","--manifest",
    dest='manifest',
    help = "Path to the manifest file",
    type = "character"
  );

  parser$add_argument(
    "-p","--platform",
    dest='platform',
    help = "Methylation array platform. Default: EPIC",
    type = "character",
    default='EPIC'
  );



  # Parse arguments
  opt <- parser$parse_args();

  for(n in names(opt)){
    if(is.null(opt[[n]])){      
      parser$print_help()
      cli::cli_alert_danger("Couldn't parse correctly the provided arguments.")
      stop()
    }
  }

  config = list()
  config[['workDir']] = opt$workDir
  config[['runName']] = opt$runName
  config[['genome']] = opt$genome
  config[['platform']] = opt$platform
  config[['manifest']] = opt$manifest
 
  return(config)
  
}



loadSourceFiles <- function(config){    

    #sourceFilesDir = normalizePath(glue::glue("{config$scriptsDir}/../../"))

    #platform = config$platform
    #genome = config$genome
        
    #manifestFile = glue::glue("{sourceFilesDir}/data/{platform}.{genome}.manifest.tsv")
    manifestFile = config$manifest
    logger::log_info("Trying to load {manifestFile}")  

    if(!file.exists(manifestFile)){
        logger::log_error("Please make sure that {manifestFile} is accessible.")
        stop()
    }

    manifest = readr::read_tsv(manifestFile, col_names=TRUE, show_col_types = FALSE) %>% as.data.frame()

    rownames(manifest) <- seq(length=nrow(manifest)) # restart row numbering    

    logger::log_info("loadSourceFiles [DONE!]")
    return(manifest)
}


prepareFolders <- function(config){

    logger::log_info("Preparing folders ...")

    dataDir = glue::glue("{config$workDir}")
    if(!dir.exists(dataDir)){
        logger::log_info("Couldn't find the data directory at {dataDir}")
        stop()
    }


    outputDir = normalizePath(".") #glue::glue("{config$workDir}/output")

    #fs::dir_create(outputDir) # This will work even if the directory exists
    #logger::log_info("Created {outputDir}")    

    # Read methyl array sheet
    logger::log_info("Reading sample sheet containing metadata ...")
    all_targets <- read.metharray.sheet(dataDir) 

    # make specific folder for this QC run
    name = config$runName
    QC_run_folder = glue::glue("{outputDir}/{name}_preprocessIllumina/")
    fs::dir_create(QC_run_folder)
    logger::log_info("Created {QC_run_folder}")

    res = list(dataDir = dataDir, 
            outputDir= QC_run_folder, 
            all_targets=all_targets
            )


    config2 = list(metharrayQCNorm = config)
    fout = glue::glue("{QC_run_folder}/config.yaml")
    yaml::write_yaml(config2,fout)

    logger::log_info("Folder preparation [DONE!]")
    return(res)
}


LoadIntoMinfi <- function(config, all_targets){

    logger::log_info("Reading methylation experiment data ...")
    RGSet.all <- read.metharray.exp(targets = all_targets, force=TRUE) # get RGChannelSet object
    logger::log_info("Methylation experiment loaded.")

    rawData_dir = glue("{config$outputDir}/raw_data")
    fs::dir_create(rawData_dir)
    
    fout=glue("{config$outputDir}/raw_data/RGSet.all.RDS") %>% normalizePath()
    logger::log_info("Saving minfi object to {fout}")
    readr::write_rds(RGSet.all,fout)

    annotation(RGSet.all) # check annotation

    phenoData.all <- pData(RGSet.all) # get the phenotype data
    phenoData.all <- as.data.frame(phenoData.all)

    fout = glue("{config$outputDir}/phenoData.all.RDS") %>% normalizePath()
    logger::log_info("Saving minfi object to {fout} ...")
    readr::write_rds(phenoData.all, fout)
    
    logger::log_info("Reading Minfi [DONE!]")
    return(RGSet.all)
}


# Original code at: https://github.com/hansenlab/minfi/blob/master/R/minfiQC.R
my.plotQC <- function(qc, badSampleCutoff = 10.5) {
    meds <- (qc$mMed + qc$uMed)/2
    whichBad <- which((meds < badSampleCutoff))
    plot(qc$mMed, qc$uMed,
         xlim = c(8,14), ylim = c(8,14), xaxt = "n", yaxt = "n",
         xlab = "Meth median intensity (log2)",
         ylab = "Unmeth median intensity (log2)",
         pch=19,
         col = ifelse(1:nrow(qc) %in% whichBad, "red", "black"))
    axis(side = 1, at = c(9,11,13))
    axis(side = 2, at = c(9,11,13))
    ## abline(h = badSampleCutoff, lty = 2)
    ## abline(v = badSampleCutoff, lty = 2)
    abline(badSampleCutoff * 2 , -1, lty = 2)
    if (length(whichBad) > 0) {
        text(qc$mMed[whichBad], qc$uMed[whichBad] - 0.25,
             labels = whichBad, col = "red")
    }
    legend("topleft", legend = c("good", "bad, with sample index"), pch = 19,
           col = c("black", "red"), bty = "n")

    qc$status = "Good"
    qc$status[whichBad] = "Bad"
    invisible(qc)
}



preNormalizationAndQC <- function(config,RGSet.all){
 
 logger::log_info("Creating un-normalized MethylSet ...")
  MSet.raw.all <- preprocessRaw(RGSet.all) # get raw methyl set
  
  logger::log_info("Estimating sample-specific quality control ...")
  qc.raw.all <- getQC(MSet.raw.all) # get quality control - methylated median and unmethylated median
  
  # get quality control plot
  fout= glue::glue("{config$outputDir}/raw_data/plotQC_raw.pdf") %>% normalizePath()
  pdf(fout)
  QC_res = my.plotQC(qc.raw.all) 
  dev.off()
  
  sampleID = pData(MSet.raw.all)$Sample_ID
  QC_res = as.data.frame(QC_res)
  rownames(QC_res) = sampleID
  QC_res = rownames_to_column(QC_res,var="Sample_ID")  
  fout = glue::glue("{config$outputDir}/raw_data/plotQC_raw_goodbadSamples.tsv")
  readr::write_tsv(QC_res, fout)

  # NOTE: should have step in here that removes failed samples; however, currently not sure how to tell which samples failed except based off of the picture with the indexes
  phenoData.all <- pData(RGSet.all) # get the phenotype data
  phenoData.all <- as.data.frame(phenoData.all)

  # plot density plots
  fout= glue::glue("{config$outputDir}/raw_data/density_beta_raw.pdf") %>% normalizePath()
  pdf(fout)
  densityPlot(MSet.raw.all, sampGroups = phenoData.all$Sample_Group)
  dev.off()

  # generate quality control report
  fout = glue::glue("{config$outputDir}/raw_data/qcReport.pdf") %>% normalizePath()
  logger::log_info("Generating QC report at: {fout}")
  qcReport(RGSet.all, pdf= fout)

  logger::log_info("Pre-normalization QC [DONE!]")  
}


sexCheck <- function(config, RGSet.all){

  logger::log_info("Mapping probs to genome ...")
  GRSet <- mapToGenome(RGSet.all) # get genomic ratio set 

  fin = glue::glue("{config$outputDir}/sex_prediction") %>% normalizePath()
  fs::dir_create(fin)
  
  logger::log_info("Sex prediction ...")
  predictedSex <- getSex(GRSet, cutoff = -2) # get the predicted sexes
  
  GRSet <- addSex(GRSet, sex = predictedSex) # add the sex prediction info to the genomic ratio set
  
  fout = glue::glue("{config$outputDir}/sex_prediction/sexplot.pdf") %>% normalizePath()
  pdf(fout) # plot the predicted sex
  plotSex(GRSet)
  dev.off()

  phenoData.all <- pData(RGSet.all) # get the phenotype data
  phenoData.all <- as.data.frame(phenoData.all)
  
  # add predicted sex information to phenoData and check to see if reported sex is equal to predicted sex
  phenoData.all$Predicted_Sex <- as.factor(predictedSex$predictedSex)
  phenoData.all$Predicted_Sex <- factor(phenoData.all$Predicted_Sex, levels = c("F", "M"),labels = c("Female","Male")) # change the factor names  
  phenoData.all$Sex_Concordance <- ifelse(as.character(phenoData.all$Predicted_Sex) == as.character(phenoData.all$Reported_Sex), TRUE, FALSE)
  

  # save sex predictions to text file
  fout = glue::glue("{config$outputDir}/sex_prediction/sexPredictionsComparison.txt")
  
  readr::write_tsv(phenoData.all[,c("Sample_Name", "Sample_ID", "Reported_Sex", "Predicted_Sex", "Sex_Concordance")], 
                   file = fout)
  
  # save sex mismatch samples to text file  
  discordant.sex <- phenoData.all %>% filter(Sex_Concordance == FALSE)
  fout= glue::glue("{config$outputDir}/sex_prediction/sexPredictionsComparison.discordant.txt")
  readr::write_tsv(discordant.sex, fout)

  logger::log_info("Sex check [DONE!]")
  return(phenoData.all)  
 
}



Normalization <- function(config, RGSet.all){

  logger::log_info("Normalizing data ...")

    
  # background subtraction and control normalization
  MSet.norm.all <- preprocessIllumina(RGSet.all, bg.correct = TRUE, normalize = "controls") 
  
  n_samples <- ncol(MSet.norm.all)
  n_probes <- nrow(MSet.norm.all)
  
  # A detection p-value is returned for every genomic position in every sample. 
  logger::log_info("Identification of failed positions")  
  detect.levels.all <- detectionP(RGSet.all, type="m+u") # The m+u method compares the total DNA signal (Methylated + Unmethylated) for each position to the background signal level. 
  
  # designate failed samples 
  failed <- detect.levels.all > 0.01 
  fraction_failed <- colMeans(failed) # fraction of failed positions per sample (# true/# total)
  fraction_failed_50 <- sum(rowMeans(failed) > 0.5) # how many positions failed in >50% of samples?

  cli::cli_alert_info("{.val {fraction_failed_50}} out of {.val {n_probes}} probes failed in >50% of the samples")

  din = glue::glue("{config$outputDir}/normalized_data")  %>% normalizePath()
  fs::dir_create(din)  
  
  
  # plot fraction of failed positions per sample
  glue::glue("Generating failed positions plot ...")
  fout = glue::glue("{din}/Fraction_Failed_Positions_PerSample.pdf")
  ttl = glue::glue("{config$platform} Fraction of Failed Positions per Sample")
  pdf(fout)
  barplot(sort(fraction_failed), main=ttl, ylab="Fraction Failed Positions")
  dev.off()
  
  # how many samples have significant p values at each site?
  detect.sd.all <- rowSums(detect.levels.all < 0.01) # sums the number of samples for which that probe[row] has a "significant" p value
  # the number associated with each probe is equivalent to the number of samples with significant p values at that site
  cli::cli_alert_info("The average number of samples with significant p-values across all probes is {.val {mean(detect.sd.all)}} out of {.val {n_samples}}.")  
  
  # which sites have >90% success rate of samples?
  detect.ind.sd.all <- which(detect.sd.all > (0.9*n_samples)) # determine which rows have > 90% of the total number of samples with "significant" p-values
  
  cli::cli_alert_info("{.val {length(detect.ind.sd.all)}} out of {.val {n_probes}} probes have >90% of the total number of samples with significant p-values")  

  detect.ind.sd.all <- as.matrix(detect.ind.sd.all) # convert to matrix
  
  # determine the overlap between the normalized methyl set and the desired sites
  overlap.all <- intersect(rownames(MSet.norm.all), rownames(detect.ind.sd.all))
  # length(overlap.all) == length(detect.ind.sd.all) # should be same length as detect.ind.sd.all
  
  target.MSet.norm.all <- MSet.norm.all[overlap.all,] # create a new methyl set containing only the desired sites (>90% success rate)
  
  number.sites.dropped <- nrow(MSet.norm.all) - nrow(target.MSet.norm.all) # check the number of probes dropped
  
  cli::cli_alert_info("Of {.val {nrow(MSet.norm.all)}} total probes, taking {.val {nrow(target.MSet.norm.all)}} probes with significant p-values in >90% of samples")  

  logger::log_info("Converting methylation signals to Beta and M-values ...")
  target.ratioSet.all <- ratioConvert(target.MSet.norm.all, what = "both", keepCN = TRUE) # convert to a ratio set

  logger::log_info("Mapping probs to genome ...")
  GSet.all <- mapToGenome(target.ratioSet.all) # map the target ratio set to the genome

  logger::log_info("Removing loci with SNPs ...")
  GSet.all <- dropLociWithSnps(GSet.all) # remove loci with snps
  
  fout = glue::glue("{config$outputDir}/normalized_data/GSet.all.RDS")
  logger::log_info("Saving  GSet.all object to : {fout}")  
  
  return(target.MSet.norm.all)
}


saveBeta <- function(MSet, phenoData, manifest, outputDir, filename) {
    
    #################################################################################################################
    # get beta values
    GSet <- mapToGenome(MSet) # convert to genomic methyl set
    GSet <- dropLociWithSnps(GSet) # remove loci with snps
    beta <- getBeta(GSet) # get beta values
    # saveRDS(beta, file = paste(outputDir, "normalized_data/beta.", filename, ".RDS", sep = "")) # save data to RDS
    
    #################################################################################################################
    # plot a density plot of the beta values
    fout = glue::glue("{outputDir}/normalized_data/beta.density.bysample.{filename}.pdf")
    pdf(fout)
    densityPlot(beta, sampGroups=phenoData$Sample_Name)
    dev.off()
    
    fout = glue::glue("{outputDir}/normalized_data/beta.density.bygroup.{filename}.pdf")
    pdf(fout)
    densityPlot(beta, sampGroups=phenoData$Sample_Group)
    dev.off()
    
    #################################################################################################################
    # annotate beta values
    # add probeID info as first column
    beta = as.data.frame(beta)
    beta = rownames_to_column(beta, var='probeID')
    #beta <- as.data.frame(cbind(rownames(beta), beta))     
    
    filtered.probes <- beta$probeID # filter the EPIC manifest data to only contain the filtered probes
    
    # filter the manifest data to only contain the filtered probes
    filtered.probes.manifest <- manifest[manifest$probeID %in% filtered.probes,] 
    filtered.probes.manifest <- filtered.probes.manifest %>% select(5, 1:3)
    
    # merge the objects at the probe IDs
    beta <- merge(filtered.probes.manifest, beta, by.x = "probeID", by.y = "probeID", all.x = FALSE, all.y = FALSE)

    # save beta values as txt file    
    fout= glue::glue("{outputDir}/normalized_data/{filename}.beta.txt")
    cli::cli_alert_info("{symbol$bullet} Saving beta values to {.val {fout}}")
    readr::write_tsv(beta, fout)    
}



seprateAutosomes <- function(configs,RGSet.all, target.MSet.norm.all, manifest, phenoData.all){

    logger::log_info("Separating autosomes from sex chromosomes ...")

    # ! = not; keep everything but sex chromosomes
    autosomes <- !(featureNames(target.MSet.norm.all) %in% manifest$probeID[manifest$CpG_chrm %in% c("chrX", "chrY")]) 
    sexchr <- (featureNames(target.MSet.norm.all) %in% manifest$probeID[manifest$CpG_chrm %in% c("chrX", "chrY")])

    MSet.norm.clean.auto <- target.MSet.norm.all[autosomes,] # autosomes for males and females

    phenoData.female <- phenoData.all[phenoData.all$Predicted_Sex == "Female",]
    females <- row.names(phenoData.female) 
    MSet.norm.clean.females.sexchr <- target.MSet.norm.all[sexchr,females] 

    phenoData.male <- phenoData.all[phenoData.all$Predicted_Sex == "Male",]
    males <- row.names(phenoData.male)
    MSet.norm.clean.males.sexchr <- target.MSet.norm.all[sexchr,males] 

    logger::log_info("Saving autosmes ...")
    saveBeta(MSet.norm.clean.auto, phenoData.all,manifest, config$outputDir, "autosomes")
    logger::log_info("Saving female chromosome ...")
    saveBeta(MSet.norm.clean.females.sexchr, phenoData.female, manifest, config$outputDir, "female.sexchr")
    logger::log_info("Saving male chromosome ...")
    saveBeta(MSet.norm.clean.males.sexchr, phenoData.male, manifest, config$outputDir, "male.sexchr")

    logger::log_info("Autosome separation [DONE!]")
}



#### 
# Helper functions for annotation

parseAnnoParameters <- function(){


  # Parse paramters
  parser <- argparse::ArgumentParser(
    description = "Main script to generated the initial QC and normalization."
  ); 

  parser$add_argument(
    "-o", "--outputDir",
    dest='outputDir',
    help = "Path to the output directory of the run",
    type = "character"
  );

  #parser$add_argument(
  #  "-n", "--runName",
  #  dest='runName',
  #  help = "Output folder name",
  #  type = "character"
  #);

  parser$add_argument(
    "-p","--platform",
    dest='platform',
    help = "Methylation array platform. Default: EPIC",
    type = "character",
    default='EPIC'
  );


  parser$add_argument(
    "-d", "--dmr",
    dest='dmr',
    help = "DMR file to annotate",
    type = "character"
  );

  parser$add_argument(
    "-g", "--genome",
    dest='genome',
    help = "Genome build to use",
    type = "character"
  );

  # Parse arguments
  opt <- parser$parse_args();

  for(n in names(opt)){
    if(is.null(opt[[n]])){      
      parser$print_help()
      cli::cli_alert_danger("Couldn't parse correctly the provided arguments.")
      stop()
    }
  }

  if(!file.exists(opt$dmr)){
      logger::log_error("Please make sure the path to the DMR file is correct.")
      cli::cli_alert_danger("Couldn't read: {.val {opt$dmr}}")
  }

  config = list()
  config[['outputDir']] = opt$outputDir
  #config[['runName']] = opt$runName
  config[['platform']] = opt$platform
  config[['genome']] = opt$genome
  config[['DMR']] = opt$dmr


  return(config)
}





runHomerAnnotation <- function(config){

    tmp_dir = tempdir(check = FALSE)    
    fin = config$DMR
    logger::log_info("Creating HOMER tmp directory at {tmp_dir}")
    fs::dir_create(tmp_dir)



    logger::log_info("Converting {config$DMR} to HOMER bed format")
    sig.df = readr::read_tsv(config$DMR, show_col_types = FALSE)
    sig.df$DMR_size = sig.df$End_DMR - sig.df$Start_DMR


    anno_dir = glue::glue("{config$outputDir}/annotations")
    fs::dir_create(anno_dir)

    fin = gsub(".txt.sig",".sig.bed",config$DMR) %>% basename()
    fin = glue::glue("{tmp_dir}/{fin}")
    readr::write_tsv(sig.df[,1:3], fin, col_names=FALSE)

    fout = gsub(".sig.bed","sig.anno.tsv",fin)
    fout_stat = gsub(".anno.tsv",".anno.stats.tsv",fout)

    annoFile = glue::glue("{anno_dir}/{basename(fout)}")
    statOut = glue::glue("{anno_dir}/{basename(fout_stat)}")



    logger::log_info("Running HOMER annotation ...")
    cmd = glue::glue("module load homer/4.10 && annotatePeaks.pl  {fin} {config$genome}> {annoFile}")
    logger::log_info("Launching HOMER with the following command: {cmd}")
    system(cmd)
    

    # Clean a little bit the results
    logger::log_info("Cleaning HOMER results ...")
    homer_res = readr::read_tsv(annoFile, show_col_types = FALSE)

    colnames(homer_res)[1] = "ID"
    homer_res = homer_res %>% arrange(ID)
    homer_res$Annotation = gsub("\\s+\\(.*","",homer_res$Annotation)
    
    # make the CpG annotation more visible
    pos = grep("CpG",homer_res[['Detailed Annotation']])
    homer_res$WithinCpGIslan = "NO"
    homer_res$WithinCpGIslan[pos] = "YES"

    # make the Repeats annotation more visible
    pos_SINE = grep("SINE",homer_res[['Detailed Annotation']])
    pos_LINE = grep("LINE",homer_res[['Detailed Annotation']])
    pos_RC = grep("RC",homer_res[['Detailed Annotation']])
    pos_LTR = grep("LTR",homer_res[['Detailed Annotation']])
    pos_LC = grep("Low_complexity",homer_res[['Detailed Annotation']])
    pos_sr = grep("Simple_repeat",homer_res[['Detailed Annotation']])
    pos_sat = grep("Satellite",homer_res[['Detailed Annotation']])


    pos = c(pos_SINE,pos_LINE, pos_RC,pos_LTR,pos_LC,pos_sr, pos_sat) %>% unique()
    homer_res$RepeatRegion = "."
    homer_res$Detailed_RepeatRegion = "."

    Repeat_info = stringr::str_match(homer_res[['Detailed Annotation']][pos],pattern = "(.+)\\|(.+)\\|(.+)")
    homer_res$Detailed_RepeatRegion[pos] = homer_res[['Detailed Annotation']][pos]
    homer_res$RepeatRegion[pos] = Repeat_info[,4]


    homer_res = cbind(sig.df,homer_res[,-c(1:7)])
    
    readr::write_tsv(homer_res, file=annoFile)

    config$anno_tmp = annoFile
    config$HOMER_stats = annoFile
    config$tmp_dir =  tmp_dir
    logger::log_info("HOMER annotation [DONE!]")
    return(config)
}


runCTCFAnnotation <- function(config){

    annoFile = glue::glue("{config$scriptsDir}/../annotations/{config$genome}/CTCF/{config$platform}.{config$genome}.CTCF.overlap.tsv")

    if(!file.exists(annoFile)){
        logger::log_error("Couldn't find CTCF annotation file.")
        cli::cli_abort("{annoFile}")
    }

    # double check the homer annotation results are still there
    if(!file.exists(config$anno_tmp)){
        logger::log_error("Couldn't find HOMER annotation file.")
        cli::cli_abort("{config$anno_tmp}")
    }

    logger::log_info("CTCF annotation ...")    
    sig.dmr.df = readr::read_tsv(config$anno_tmp, show_col_types = FALSE)
    old_cnames = colnames(sig.dmr.df)
    colnames(sig.dmr.df)[1:3] = c("seqnames","start","end")
    sig.dmr.gr = GRanges(sig.dmr.df)

    colnames(sig.dmr.df)[1:3] = old_cnames[1:3]

    # here readr::read_tsv doesn't work somehow
    ctcf_anno = read.table(annoFile, sep="\t",row.names=1)
    colnames(ctcf_anno)[1:3] = c("seqnames","start","end")

    ctcf_anno.gr = GRanges(ctcf_anno)

    ovp = findOverlaps(sig.dmr.gr, ctcf_anno.gr)

    dmr.hits = table(queryHits(ovp))

    sig.dmr.df$nOverlapping_CTCFprobs = 0
    sig.dmr.df$nOverlapping_CTCFprobs[as.numeric(names(dmr.hits))] = dmr.hits

    logger::log_info("Updating annoration ...")  
    readr::write_tsv(sig.dmr.df, file=config$anno_tmp)

    logger::log_info("CTCF annotation [DONE!]")
    return(config)
}


runImprintingAnnotation <- function(config){

  annoFile = glue::glue("{config$scriptsDir}/../annotations/{config$genome}/imprinting/{config$platform}.{config$genome}.imprinting.tsv")

  if(!file.exists(annoFile)){
      logger::log_error("Couldn't find imprinting annotation file.")
      cli::cli_abort("{annoFile}")
  }

  # double check the homer annotation results are still there
  if(!file.exists(config$anno_tmp)){
      logger::log_error("Couldn't find HOMER annotation file.")
      cli::cli_abort("{config$anno_tmp}")
  }

  logger::log_info("Imprinting annotation ...")    
  sig.dmr.df = readr::read_tsv(config$anno_tmp, show_col_types = FALSE)
  colnames(sig.dmr.df)[1:3] = c("seqnames","start","end")
  sig.dmr.gr = GRanges(sig.dmr.df)

  # here readr::read_tsv doesn't work somehow
  imprinting_anno = readr::read_tsv(annoFile)
  colnames(imprinting_anno)[1:3] = c("seqnames","start","end")

  imprinting_anno.gr = GRanges(imprinting_anno)

  ovp = findOverlaps(sig.dmr.gr, imprinting_anno.gr)

  dmr.hits = table(queryHits(ovp))

  sig.dmr.df$nOverlapping_Imprintingprobs = 0
  sig.dmr.df$nOverlapping_Imprintingprobs[as.numeric(names(dmr.hits))] = dmr.hits

  dmr.patmat = as.data.frame(ovp)
  colnames(dmr.patmat)= c("DMR_ID","Prob_ID")

  dmr.patmat$PaternalOrMaternal = imprinting_anno.gr$PaternalOrMaternal[dmr.patmat$Prob_ID]
  

  dmr.patmat_stats = split(dmr.patmat, dmr.patmat$DMR_ID) %>%
      imap(.f= function(df, DMR_id){

        res= tibble(DMR_ID = DMR_id, 
                    PaternalOrMaternal = paste(sort(unique(df$PaternalOrMaternal)),collapse=",")
                    )
        return(res)
      }) %>%
      do.call(what="rbind")


  sig.dmr.df$PaternalOrMaternal = "."
  sig.dmr.df$PaternalOrMaternal[as.numeric(dmr.patmat_stats$DMR_ID)] = dmr.patmat_stats$PaternalOrMaternal

  logger::log_info("Updating annoration ...")  
  readr::write_tsv(sig.dmr.df, file=config$anno_tmp)

  logger::log_info("Imprinting annotation [DONE!]")
  return(config)
}


GRangesToString <- function (grange, sep = c(":", "-")) 
{
    regions <- paste0(as.character(x = seqnames(x = grange)), 
        sep[[1]], start(x = grange), sep[[2]], end(x = grange))
    return(regions)
}

run450KAnno <- function(cofig){

  annoFile = glue::glue("{config$scriptsDir}/../annotations/{config$genome}/450k/450k_population_DMRs_autosomes.xlsx")

  if(!file.exists(annoFile)){
      logger::log_error("Couldn't find imprinting annotation file.")
      cli::cli_abort("{annoFile}")
  }

  # double check the homer annotation results are still there
  if(!file.exists(config$anno_tmp)){
      logger::log_error("Couldn't find annotation file.")
      cli::cli_abort("{config$anno_tmp}")
  }

  logger::log_info("Population frequency annotation ...")    
  sig.dmr.df = readr::read_tsv(config$anno_tmp, show_col_types = FALSE)
  org_names = colnames(sig.dmr.df)[1:3]
  colnames(sig.dmr.df)[1:3] = c("seqnames","start","end")  
  sig.dmr.gr = GRanges(sig.dmr.df)

  colnames(sig.dmr.df)[1:3] = org_names
  sig.dmr.df[["Frequency of DMR per 10,000 (95% CI)"]] = '.'
  sig.dmr.df[["Corresponding 450k population DMR"]] = '.'

  # here readr::read_tsv doesn't work somehow
  logger::log_info("Loading 450k population data ...")
  population_anno = readxl::read_xlsx(annoFile, skip=2)
  colnames(population_anno)[1:3] = c("seqnames","start","end")

  population_anno.gr = GRanges(population_anno)

  sig.dmr.gr$direction = gsub("MAX","HYPER",sig.dmr.gr$direction)
  sig.dmr.gr$direction = gsub("MIN","HYPO",sig.dmr.gr$direction)

  sig.dmr.gr$rowID = 1:length(sig.dmr.gr)

  sig.dmr.gr_lst = split(sig.dmr.gr, sig.dmr.gr$direction)


  ## simpligy the dirrection of the 450k samples
  dir = sapply(mcols(population_anno.gr)[["Direction of methylation change"]], 
                function(x) {
                  splt = strsplit(x,split=',') %>% unlist()
                  splt = gsub("\\s+","",splt) %>% unique()
                  if(length(splt)>1){
                    return("Both")
                  }
                  splt = ifelse(splt =="Gain","HYPER","HYPO")
                  return(splt)
                })

  population_anno.gr$Direction = dir

  population_anno.gr_lst = split(population_anno.gr, population_anno.gr$Direction)

  ## Overlap only DMR with the same directoin
  for(direction in names(sig.dmr.gr_lst)){
    
    ovp = findOverlaps(sig.dmr.gr_lst[[direction]], population_anno.gr_lst[[direction]])

    if(length(ovp)>0){      
      dmr.hits = table(queryHits(ovp))

      unique_hit_dmr = names(dmr.hits[dmr.hits == 1])
      multi_hit_dmr = names(dmr.hits[dmr.hits>1])      
      
      unique_ovp = subset(as.data.frame(ovp), queryHits %in% unique_hit_dmr)
      
      DMR_id = sig.dmr.gr_lst[[direction]]$rowID[unique_ovp$queryHits ]

      sig.dmr.df[["Frequency of DMR per 10,000 (95% CI)"]][DMR_id] = mcols(population_anno.gr_lst[[direction]])[["Frequency of DMR per 10,000 (95% CI)"]][unique_ovp$subjectHits]      
      sig.dmr.df[["Corresponding 450k population DMR"]][DMR_id] = GRangesToString(population_anno.gr_lst[[direction]][unique_ovp$subjectHits])

      # fix multi-hits DMR if any 
      if(length(multi_hit_dmr)>0){
        logger::log_info("{length(multi_hit_dmr)} {direction} DMR(s) overlapped with multiple 450k DMRs.")  

        for(dmrid in multi_hit_dmr){
          gr = subsetByOverlaps(population_anno.gr_lst[[direction]], sig.dmr.gr_lst[[direction]][as.numeric(dmrid)])       
          intervals_df = stringr::str_match(mcols(gr)[["Frequency of DMR per 10,000 (95% CI)"]],pattern = "(.*)\\s*-\\s*(.*)")
          max_interval = apply(intervals_df[,2:3],2,max)

          DMR_ID = sig.dmr.gr_lst[[direction]][as.numeric(dmrid)]$rowID
          # We mark these regions with an asterics 
          sig.dmr.df[["Frequency of DMR per 10,000 (95% CI)"]][DMR_ID] = glue::glue("{max_interval[1]} - {max_interval[2]}*")

          overlapping_dmrs = GRangesToString(gr) 
          overlapping_dmrs = paste(overlapping_dmrs, collapse = ",")
          sig.dmr.df[["Corresponding 450k population DMR"]][DMR_ID] = overlapping_dmrs
        }
      }      
    }    
  }  

  logger::log_info("Updating annoration ...")  
  readr::write_tsv(sig.dmr.df, file=config$anno_tmp)

  logger::log_info("Population frequency annotation [DONE!]")
  return(config)

}

saveAnnoAsXlsx <- function(config){

  # double check the homer annotation results are still there
    if(!file.exists(config$anno_tmp)){
        logger::log_error("Couldn't find HOMER annotation file.")
        cli::cli_abort("{config$anno_tmp}")
    }

  
  sig.dmr.df = readr::read_tsv(config$anno_tmp, show_col_types = FALSE)

  fout = gsub(".tsv",".xlsx",config$anno_tmp)
  logger::log_info(glue("Saving results to {fout} ..."))
  wb <- createWorkbook()
  addWorksheet(wb = wb,sheetName = "Marker genes")

  writeDataTable(wb,sheet = 1,
                 x = sig.dmr.df,
                 tableStyle ="TableStyleLight15",
                 bandedCols = F,
                 bandedRows = F)


  hs1 <- createStyle(fgFill = "#F2F2F2",fontColour = "black",textRotation = 60,textDecoration = "Bold")
  addStyle(wb, 1, style = hs1, rows = 1, cols=1:ncol(sig.dmr.df))

  fout = gsub(".tsv",".xlsx",config$anno_tmp)
  saveWorkbook(wb,file = fout, overwrite = TRUE)

}
