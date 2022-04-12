
###############################################################################
# Load dependencies and helper scripts

## Check if funr is installed (it is essential to get the source file path).
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
cli::cli_blockquote("DMR annotation")

cli::cli_h1("Loading and installing missing packages")
checkPackages()



###############################################################################
# Parse parameters
cli::cli_h1("Parse parameters and folder prep")
config  = parseAnnoParameters()
config$scriptsDir = rootDir

#readr::write_rds(config, "config.rds")


###############################################################################
# Run HOMER annotation
cli::cli_h1("Running Annotation")
cli::cli_h3("HOMER annotation ..")
config = runHomerAnnotation(config)

cli::cli_h3("CTCF annotation ...")
config = runCTCFAnnotation(config)


cli::cli_alert_info("Imprinting genes annotation ...")
config = runImprintingAnnotation(config)


cli::cli_alert_info("Population Frequency annotation ...")
config = run450KAnno(config)

readr::write_rds(config, file="config.rds")

cli::cli_alert_info("Exporting annotation to xlsx format ...")
saveAnnoAsXlsx(config)

