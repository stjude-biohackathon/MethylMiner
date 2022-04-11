# **SJHack2022_Project2**
A methylation array analysis pipeline tailored for discovering rare methylation events with interactive data visualization.

# **Team Members**

    - Christy LaFlamme
    - Pandurang Kolekar
    - Nadhir Djekidel
    - Wojciech Rosikiewicz

Methylation array analysis pipeline consists of two major components

1. Processing of data using NextFlow pipeline
2. Interactive exploration using Jupyter Dash server

# **Installation**
 
## **1. NextFlow Pipeline**

Current implementation of the pipeline is designed to be executed on St. Jude HPC Research Cluster.

After logging into interactive node, please load the following NextFlow module before running the pipeline

```
module load nextflow/21.04.1
```

### **Pipeline input parameters**

* **workdir**    - Input directory with *.idat files and single *.csv file with description of the samples
* **runName**    - Name for the execution run of the pipeline
* **outdir**     - Output directory for the pipeline to save the output files/folders
* **genome**     - genome for QC/DMR analysis (currently only support hg19)
* **platform**   - Methylation array platform (currently only support EPIC)
* **windowSize** - Window size (bp) for DMR analysis (default: 1000)
* **qCutMin**    - Minimum quantile cut-off for DMR analysis (default: 0.25)
* **qCutMax**    - Maximum quantile cut-off for DMR analysis (default: 0.75)
* **email**      - Email address for notification of the pipeline status

### **.csv file datasheet REQUIRED columns**

* **Sample_Name**       - Column containing the sample names
* **Sentrix_ID**        - Should be first portion of unique identifier consisting of 12 numbers (e.g. "204776850065") - Be careful! Microsoft Excel likes to convert these to scientific notation
* **Sentrix_Position**  - Should be row column number (e.g. R01C01)
* **Sample_ID**         - Should be the Sentrix_ID and the Sentrix Position corresponding to the full unique identifier (e.g. "204776850065_R01C01") 
* **Reported_Sex**      - Options: Male, Female, Unknown (required if sex check for QC is desired)
* **Sample_Group**      - Option: Case, Control (required if post-analysis filtering is desired)

Columns must be named exactly as indicated above. Also, .csv file must end on a new line (return character).

### **.csv file datasheet RECOMMENDED columns**

* **Sample_Type**       - Column containing the sample tissue source (e.g. Blood, Saliva, Fibroblast, etc.)
* **Sentrix_Well**      - Well number corresponding to the sample plate submission
* **Sample_Plate**      - Plate name that sample was run on (good annotation to check for batch effects between plates)
* **Date**              - Date of sample run

### **Pipeline output**

* **Directory-level output organization**

```
<outdir>
└── <runName>
    └── <runName>_preprocessIllumina
        ├── normalized_data
        │   ├── betaValues
        │   ├── bigWig
        │   │   ├── autosomes
        │   │   ├── female
        │   │   └── male
        │   ├── dmr_<windowSize>_<qCutMin>_<qCutMax>
        │   └── qcReports
        ├── raw_data
        └── sex_prediction
```

* **Directory/File-level output organization**

```
<outdir>
└── <runName>
    └── <runName>_preprocessIllumina
        ├── normalized_data
        │   ├── betaValues
        │   │   ├── autosomes.beta.txt.sorted
        │   │   ├── female.sexchr.beta.txt.sorted
        │   │   └── male.sexchr.beta.txt.sorted
        │   ├── bigWig
        │   │   ├── autosomes
        │   │   │   ├── <sample_1>.bed
        │   │   │   ├── <sample_1>.bw
        │   │   │   ├── <sample_2>.bed
        │   │   │   ├── <sample_2>.bw
        │   │   │   ├── <sample_n>.bed
        │   │   │   └── <sample_n>.bw
        │   │   ├── female
        │   │   │   ├── <sample_1>.bed
        │   │   │   ├── <sample_1>.bw
        │   │   │   ├── <sample_2>.bed
        │   │   │   ├── <sample_2>.bw
        │   │   │   ├── <sample_n>.bed
        │   │   │   └── <sample_n>.bw
        │   │   └── male
        │   │       ├── <sample_1>.bed
        │   │       ├── <sample_1>.bw
        │   │       ├── <sample_2>.bed
        │   │       ├── <sample_2>.bw
        │   │       ├── <sample_n>.bed
        │   │       └── <sample_n>.bw
        │   ├── dmr_<windowSize>_<qCutMin>_<qCutMax>
        │   │   ├── [autosomes|sexchr].beta.txt.sorted_<windowSize>_<qCutMin>_<qCutMax>_findEpivariation.txt
        │   │   ├── [autosomes|sexchr].beta.txt.sorted_<windowSize>_<qCutMin>_<qCutMax>_findEpivariation.txt.sig
        │   │   └── [autosomes|sexchr].beta.txt.sorted_<windowSize>_<qCutMin>_<qCutMax>_findEpivariation.txt.sig_dmr.txt
        │   └── qcReports
        │       ├── beta.density.bygroup.[autosomes|sexchr].pdf
        │       └── beta.density.bysample.[autosomes|sexchr].pdf
        ├── phenoData.all.RDS
        ├── raw_data
        │   ├── density_beta_raw.pdf
        │   ├── plotQC_raw.pdf
        │   ├── qcReport.pdf
        │   └── RGSet.all.RDS
        └── sex_prediction
            ├── sexplot.pdf
            ├── sexPredictionsComparison.discordant.txt
            └── sexPredictionsComparison.txt
```

### **How to run the pipeline?**

Clone the repository on your cluster directory and navigate into pipeline directory

Pipeline consists of the following three workflows.


- **QC_WF**
    * Workflow to perform only Quality Checks (QC) on input raw data
    * Command format:
    
        ```
        nextflow run map.nf -entry QC_WF \
            --workdir <Path/to/WorkDir> \
            --runName <RunName> \
            --outdir <Path/to/OutDir> \
            --email <email address>
        ```
    
- **DMR_WF**
    * Workflow to perform Differentially Methylated Region (DMR) analysis on output generated by _QC_WF_ workflow
    * _Command format_:
    
    
        ```
        nextflow run map.nf -entry DMR_WF \
            --workdir <Path/to/WorkDir> \
            --runName <RunName> \
            --outdir <Path/to/OutDir> \
            --windowSize <num> \
            --qCutMin <qCutMin> \
            --qCutMax <qCutMin> \
            --email <email address>
        ```
        
    * **Note**
    
        - _DMR_WF_ workflow can be run multiple times on the same QC output with different window size and quantile cut-off parameters.
        - Pipeline will generate separate output directories prefixing the related parameters
    
- **QUANTILE_WF**
    * Workflow to perform QC abd DMR analysis on input raw data
    * _Command format_:
    
    
        ```
        nextflow run map.nf -entry DMR_WF \
            --workdir <Path/to/WorkDir> \
            --runName <RunName> \
            --outdir <Path/to/OutDir> \
            --windowSize <num> \
            --qCutMin <qCutMin> \
            --qCutMax <qCutMin> \
            --email <email address>
        ```

### **Future updates to the pipeline**
    - Support for different genomes (hg38, mm9, mm10)
    - Support for different methylation array platforms (450K etc)
    - Support for different normalization methods
    - Support for cell type composition QC analysis
    - Support for CNV calling using "conumee"
    - Docker/Singularity containerization for cross-platform execution with required dependencies

## Jupyter Dash 

# Setup
Directory setup (workDir)
data folder, .idat files, csv (required areas of CSV)

# Features
-Quality Control and Normalization


input:
output:

-Calling Rare Differentially Methylated Regions (DMR) 


input:
output:

-Interactive Data Visualization


input:
output:

# Usage/Tutorial
Download sample/trial data

