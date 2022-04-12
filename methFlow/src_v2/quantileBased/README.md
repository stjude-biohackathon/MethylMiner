
This is just a refectored version of the initial R scripts. The goal is to make it more configurable and flexible. 
All the parameters are stored in config file.


# QC normalization script

The QC normalization script loads the data to `minfi` and generates several QC parameters, it also normalizes the data and separates automsome and sex chromosomes.

## Command line parameters
```bash
usage: 1_metharray_QC_norm.R [-h] [-w WORKDIR] [-n RUNNAME] [-g GENOME]
                             [-m MANIFEST] [-p PLATFORM]

Main script to generated the initial QC and normalization.

optional arguments:
  -h, --help            show this help message and exit
  -w WORKDIR, --workDir WORKDIR
                        Path to the working directory
  -n RUNNAME, --runName RUNNAME
                        Output folder name
  -g GENOME, --genome GENOME
                        Genome build. Default: hg19
  -m MANIFEST, --manifest MANIFEST
                        Path to the manifest file
  -p PLATFORM, --platform PLATFORM
                        Methylation array platform. Default: EPIC
```


## Example

```bash
module load R/4.1.0

workDir="~/hackathlon/2022_methylation_array_biohackathon_challenge/"
runName="TestRun"
genome="hg19"
platform="EPIC"
manifest="/research/rgs01/home/clusterHome/mdjekide/hackathlon/2022_methylation_array_biohackathon_challenge/source_files/EPIC.hg19.manifest.tsv"

Rscript 1_metharray_QC_norm.R -w ${workDir} -n ${runName} -g ${genome} -p ${platform} -m ${m}
```

# Annotation script

The goal of this script is to annotate the identified DMR regions to different annotations we have.
We use HOMER to annotate: Genomic Regions, CpG island and repeats.

We use the plateform assocated CTCF and Imprinting files to added mark any DMRs overlapping CTCF or the ones located in imprinted regions.

DMRs are also compared to large scale DMRs of normal patients obtaine from `XXXX et al, ....` to get the population frequency. Generally a DMR will overlap one reference DMR, however, few will overlap with multiple regions. These are marked with `*` in the final annotation file.

## Command line parameters

```bash
usage: 5_annotateDMRs.R [-h] [-o OUTPUTDIR] [-p PLATFORM] [-d DMR] [-g GENOME]

Main script to generated the initial QC and normalization.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUTDIR, --outputDir OUTPUTDIR
                        Path to the output directory of the run
  -p PLATFORM, --platform PLATFORM
                        Methylation array platform. Default: EPIC
  -d DMR, --dmr DMR     DMR file to annotate
  -g GENOME, --genome GENOME
                        Genome build to use
```

## Example

```bash
module load R/4.1.0

outDir="/home/mdjekide/hackathlon/2022_methylation_array_biohackathon_challenge/output/10_test_metharray_preprocessIllumina/"
genome="hg19"
platform="EPIC"
DMR_sig="~/tmp/biohackathon_challenge/output/10_test_metharray_preprocessIllumina/normalized_data/autosomes.beta_dmr.txt.sorted.0.75.0.25.findEpivariation.txt.sig"

Rscript 5_annotateDMRs.R -o ${outDir}  -p ${platform} -d ${DMR_sig} -g ${genome} -p ${platform}
```