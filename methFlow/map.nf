#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   St. Jude BioHackathon 2022
   Team 2 : Methylation array analysis pipeline (MAP)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/stjude-biohackathon/SJHack2022_Project2
    Slack  : https://stjude.slack.com/archives/C036KD748GZ    
----------------------------------------------------------------------------------------
*/
    
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    METHYLATION PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.workdir      = "${projectDir}/test/testWorkDir"
params.runName      = 'testRun'
params.outdir       = "${projectDir}/test/testOutDir"
params.genome       = 'hg19'
params.platform     = 'EPIC'
params.windowSize   = 1000
params.qCutMin      = 0.25
params.qCutMax      = 0.75
params.email        = null

process QC_NORMALIZATION {
    publishDir "${params.outdir}/${params.runName}", mode: 'copy'
    cache 'lenient'
    //executor 'lsf'
    //memory '8 GB'
    //queue 'compbio'
    
    input:
        val wdir
        val rName
        val plfrm
    output:
        path "${rName}_preprocessIllumina/", emit: outpath
    
    """
    module load R/4.1.0
    manifest="${projectDir}/data/EPIC.hg19.manifest.tsv"
    if [ ! -f \${manifest} ]; then
        wget -P "${projectDir}/data/" "https://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/EPIC/EPIC.hg19.manifest.tsv.gz"
        gunzip "${projectDir}/data/EPIC.hg19.manifest.tsv.gz"
    fi
    script="${projectDir}/src_v2/quantileBased/1_metharray_QC_norm.R"
    cwd=\$(realpath ./)
    Rscript --vanilla \${script} -w ${wdir} -n ${rName} -m \${manifest} -p ${plfrm}
    mkdir "${rName}_preprocessIllumina/normalized_data/qcReports/"
    csvFile=\$(realpath ${wdir}/*.csv)
    mv ${rName}_preprocessIllumina/normalized_data/*.pdf ${rName}_preprocessIllumina/normalized_data/qcReports/
    cp \$csvFile ${rName}_preprocessIllumina/raw_data/
    rdsFile=\$(find ${rName}_preprocessIllumina -name RGSet.all.RDS)
    rm \$(rdsFile) || true
    """
}

process SORT_BETA {
    publishDir "${params.outdir}/${params.runName}/${dir}/normalized_data/betaValues", mode: 'copy'
    echo true
    cache 'lenient'
    //executor 'lsf'
    //memory '8 GB'
    //queue 'compbio'
    
    input:
        val out

    output:
        path "*.sorted", emit: sortBeta
    
    script:
        dir = file(out).getBaseName()
    
    """
    for FILE in "$out/normalized_data/"*.beta.txt
    do 
        head -n 1 \${FILE} > "\${FILE}.sorted" ; tail -n +2 \${FILE}  | sort -k2,2 -k3,3g -k4,4g >> "\${FILE}.sorted"
        cp "\${FILE}.sorted" ./
    done
    rm ${params.outdir}/${params.runName}/${dir}/normalized_data/*.txt
    """
}

process BED_BIGWIG {
    publishDir "${params.outdir}/${params.runName}/${params.runName}_preprocessIllumina/normalized_data/bigWig", mode: 'copy'
    publishDir "${params.outdir}/${params.runName}/${params.runName}_preprocessIllumina/assets/", mode: 'copy'
    echo true
    cache 'lenient'
    //executor 'lsf'
    //memory '8 GB'
    //queue 'compbio'
    
    input:
        path sortedBeta
    output:
        path "${dirOut}/"
    
    script:
    dirOut = file(sortedBeta).getSimpleName()
    
    """
    module load python/3.7.0
    script="${projectDir}/src/quantileBased/getBigWig.py"
    infile=\$(realpath *.sorted)
    python \${script} -betaFile \${infile} -exePath "${projectDir}/bin/bedGraphToBigWig" -chrFile "${projectDir}/data/hg19.chrom.sizes"
    
    """


}

process FIND_EPIVARIATION {
    publishDir "${params.outdir}/${params.runName}/${params.runName}_preprocessIllumina/normalized_data/dmr_${params.windowSize}_${params.qCutMin}_${params.qCutMax}", mode: 'copy'
    //echo true
    cache 'lenient'
    //executor 'lsf'
    //memory '8 GB'
    //queue 'compbio'
    
    input:
        path sortedFile
    output:
        path "*_findEpivariation.txt", emit: epiFile
        
    """
    module load perl/5.10.1
    script="${projectDir}/src/quantileBased/2_findEpivariation.pl"
    sample_end=\$((\$(head -n 1 ${sortedFile} | wc -w) + 1))
    outfile="${sortedFile}_${params.windowSize}_${params.qCutMax}_${params.qCutMin}_findEpivariation.txt"
    perl \${script} -f ${sortedFile} -c 1 -s 2 -e 3  -w ${params.windowSize} -p 3 -a 4 -b \$sample_end -1 ${params.qCutMin} -9 ${params.qCutMax} > \$outfile
    
    """
}

process FIND_DMRS {
    publishDir "${params.outdir}/${params.runName}/${params.runName}_preprocessIllumina/normalized_data/dmr_${params.windowSize}_${params.qCutMin}_${params.qCutMax}", mode: 'copy'
    //echo true
    cache 'lenient'
    //executor 'lsf'
    //memory '8 GB'
    //queue 'compbio'
    
    input:
        path epiVarFile
    output:
        path "*_findEpivariation.txt.sig", emit: epiSigFile
        path "*_dmr.txt" optional true
    
    """
    module load R/4.1.0
    script="${projectDir}/src/quantileBased/3_getDMRlist.R"
    sig_col=\$((\$(head -n 1 ${epiVarFile} | wc -w)))
    head -n 1 ${epiVarFile} > "${epiVarFile}.sig"
    cat ${epiVarFile} | awk -v NUM=\$sig_col '\$NUM == 1' >> "${epiVarFile}.sig"
    
    sigNum=\$(cat "${epiVarFile}.sig" | wc -l)
    if [[ \$sigNum > 1 ]]
    then
        Rscript \${script} "${epiVarFile}.sig" ${params.qCutMin} ${params.qCutMax}
    fi
    
    """

}


process ANNOTATE_DMR {
    publishDir "${params.outdir}/${params.runName}/${params.runName}_preprocessIllumina/normalized_data/dmr_${params.windowSize}_${params.qCutMin}_${params.qCutMax}", mode: 'copy'
    //echo true
    cache 'lenient'
    //executor 'lsf'
    //memory '8 GB'
    //queue 'compbio'
    
    input:
        path epiVarFile
    output:
        path "*_findEpivariation.txt.sig", emit: epiSigFile
        path "*_dmr.txt" optional true

    """
    module load R/4.1.0
    script="${projectDir}/src/quantileBased/3_getDMRlist.R"

    """
}


workflow QC_WF {
    wch = Channel.of(params.workdir)
    run = Channel.of(params.runName)
    plfrm= Channel.of(params.platform)
    QC_NORMALIZATION(wch, run,plfrm)
    SORT_BETA(QC_NORMALIZATION.out.outpath)
    //SORT_BETA.out.flatMap().view()
    BED_BIGWIG(SORT_BETA.out.flatten())
}

workflow DMR_WF {
    sortCh = Channel.fromPath("${params.outdir}/${params.runName}/${params.runName}_preprocessIllumina/normalized_data/betaValues/*.sorted")
    FIND_EPIVARIATION(sortCh)
    FIND_DMRS(FIND_EPIVARIATION.out.epiFile)
}

workflow QUANTILE_WF {
    wch = Channel.of(params.workdir)
    run = Channel.of(params.runName)
    QC_NORMALIZATION(wch, run)
    SORT_BETA(QC_NORMALIZATION.out.outpath)
    BED_BIGWIG(SORT_BETA.out.flatten())
    FIND_EPIVARIATION(SORT_BETA.out.flatten())
    FIND_DMRS(FIND_EPIVARIATION.out.epiFile)
}

workflow.onComplete {

    def msg = """\
    
        Methylation array analysis pipeline
        ---------------------------
        Work dir      : ${params.workdir}
        Run name      : ${params.runName}
        Genome        : ${params.genome}
        Platform      : ${params.platform}
        Window size   : ${params.windowSize}
        qCutMin       : ${params.qCutMin}
        qCutMax       : ${params.qCutMax}        
        Output dir    : ${params.outdir}
        Completed at  : ${workflow.complete}
        Duration      : ${workflow.duration}
        Success       : ${workflow.success}
        launchDir     : ${workflow.launchDir}
        exit status   : ${workflow.exitStatus}
        """
        .stripIndent()
    println msg

    sendMail(to: "${params.email}", subject: "Methylation array pipeline: ${params.runName}", body: msg)
}
