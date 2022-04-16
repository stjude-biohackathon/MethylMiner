![MethMiner](./data/assets/MethMiner.png)

# **MethMiner**
A methylation array analysis pipeline tailored for discovering rare methylation events with interactive data visualization.

A differentially methylated region (DMR) is a genomic region that has different DNA methylation status across samples. Rare DMRs have been shown to cause multiple diseases. A famous example is Fragile X sydrome where the majority of patients display a hypermethylated DMR in the 5'UTR of *FMR1* leading to reduced gene expression. The hypermethylation is driven by an underying GC-rich expansion repeat that is impossible to detect by short-read sequencing approaches. Therefore, it is very useful to "mine" these rare methylation events in order to discover novel causes of disease. 

# **Team Members**

    - Christy LaFlamme
    - Pandurang Kolekar
    - Nadhir Djekidel
    - Wojciech Rosikiewicz

Methylation array analysis pipeline consists of two major components:

1. Processing of data using NextFlow pipeline
2. Interactive exploration using Jupyter Dash server

# **Installation and usage**
 
## **1. methFlow - NextFlow Pipeline**

Current implementation of the pipeline is designed to be executed on St. Jude HPC Research Cluster.

After logging into interactive node, please load the following or latest available NextFlow module before running the pipeline

```
module load nextflow/21.10.5
```

### **Pipeline input parameters**

* **workdir**    - Input directory with *.idat files and single *.csv (description in below section) file with description of the samples
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
* **Sentrix_ID**        - Should be first portion of unique identifier consisting of 12 numbers (e.g. "204776850065")
* **Sentrix_Position**  - Should be row column number (e.g. R01C01)
* **Sample_ID**         - Should be the Sentrix_ID and the Sentrix Position corresponding to the full unique identifier (e.g. "204776850065_R01C01") 
* **Reported_Sex**      - Options: Male, Female, Unknown (required if sex check for QC is desired)
* **Sample_Group**      - Option: Case, Control (required if post-analysis filtering is desired)

Columns must be named exactly as indicated above. Also, .csv file must end on a new line (return character). And be careful! Microsoft Excel likes to convert Sentrix IDs to scientific notation.

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

Clone the repository on your cluster directory and navigate into a pipeline directory, methFlow.

The default cluster parameters in the ```nextflow.config``` file can be modified suitably. Please refer to the configuration section of [NextFlow Documentation](https://www.nextflow.io/docs/latest/index.html) for the same.

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
        - Use the same runName and outdir parameters as for _QC_WF_
        - Pipeline will generate separate output directories prefixing the related parameters
    
- **QUANTILE_WF**
    * Workflow to perform QC abd DMR analysis on input raw data
    * _Command format_:
    
    
        ```
        nextflow run map.nf -entry QUANTILE_WF \
            --workdir <Path/to/WorkDir> \
            --runName <RunName> \
            --outdir <Path/to/OutDir> \
            --windowSize <num> \
            --qCutMin <qCutMin> \
            --qCutMax <qCutMin> \
            --email <email address>
        ```

### **How to visualize and explore the results of the pipeline?**

Place the `app.py` file from MethVis subfolder inside the `<runName>_preprocessIllumina` directory, activate the environment with all dash-related packages (see installation instructions above), and run command `python app.py`.


### **Future updates to the pipeline**
    - Support for different genomes (hg38, mm9, mm10)
    - Support for different methylation array platforms (450K etc)
    - Support for different normalization methods (preprocessNoob etc)
    - Support for cell type composition QC analysis (estimateCellCounts)
    - Support for CNV calling using (conumee)
    - Docker/Singularity containerization for cross-platform execution with required dependencies

## **2. MethVis - Jupyter Dash server installation**

We recommend generating the [Anaconda](https://www.anaconda.com/products/distribution) virtual environment, to install the required packages. Assuming this is the approach taken, first create and activate the enviroment:

```
conda create -n dash
conda activate dash
```

then install all dependencies inside:

```
conda install jupyterLab numpy scipy matplotlib seaborn pandas matplotlib-venn ipykernel plotly dash dash-bio dash-core-components dash-bootstrap-components notebook jupyter-dash -c conda-forge -c anaconda -c plotly
```

### How to run the MethVis server?

```
$ cp ./MethVis/app.py /path/to/<runName>_preprocessIllumina/
$ cd /path/to/<runName>_preprocessIllumina/ 
$ conda activate dash
$ python app.py
``` 

The script will process the output files and output URL for the hosting Jupyter Dash server.

### **Future updates to the Jupyter Dash server**
    - Option to automatically run the server at the end of pipeline execution
    - More annotation tracks and interactive visualization plots

### **Test data**
Six (6) .idat files are made available as test data for the MethMiner pipeline. These samples serve as positive controls for DMR calling and sample data visualization. Three (3) samples (2 male Fragile X patients, 1 female Fragile X patient) are positive controls for a hypermethylated DMR upstream *FMR1* gene on ChrX. Three (3) samples (two affected Baratela-Scott syndrome probands, one unaffected parent) are positive controls for a hypermethylated DMR on Chr16. The metadata for these samples can be found in the same folder as the data.

1. [Fragile X NA09145](https://www.coriell.org/0/Sections/Search/Sample_Detail.aspx?Product=DNA&Ref=NA09145)
2. [Fragile X NA09237](https://www.coriell.org/0/Sections/Search/Sample_Detail.aspx?Product=DNA&Ref=NA09237)
3. [Fragile X NA07063](https://www.coriell.org/0/Sections/Search/Sample_Detail.aspx?Product=DNA&Ref=NA07063)
4. [Baratela-Scott syndrome](https://pubmed.ncbi.nlm.nih.gov/30554721/)

### **Citations**

> **_DMR calling_** : Garg P, Jadhav B, Rodriguez OL, Patel N, Martin-Trujillo A, Jain M, Metsu S, Olsen H, Paten B, Ritz B, Kooy RF, Gecz J, Sharp AJ. A Survey of Rare Epigenetic Variation in 23,116 Human Genomes Identifies Disease-Relevant Epivariations and CGG Expansions. Am J Hum Genet. 2020;107(4):654-69. Epub 2020/09/17. doi: 10.1016/j.ajhg.2020.08.019. PubMed PMID: 32937144; PMCID: PMC7536611.

> **_DMR annotations_** : Zhou W, Laird PW, Shen H. Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes. Nucleic Acids Res. 2017;45(4):e22. Epub 2016/12/08. doi: 10.1093/nar/gkw967. PubMed PMID: 27924034; PMCID: PMC5389466.
