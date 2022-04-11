#!/usr/bin/env python
# coding: utf-8

# ### Generate bed and bigWig files from sorted beta value files

import os
import sys
import argparse
import math
from glob import glob
import pandas as pd


def getSamples(betaDf):
    defCols = ['probeID', 'CpG_chrm', 'CpG_beg', 'CpG_end']
    return [col for col in betaDf.columns if col not in defCols]


def getMvalPerSample(betaDf, sampleName, mval='True'):
    defCols = ['CpG_chrm', 'CpG_beg', 'CpG_end']
    betaDf.index = betaDf.probeID
    sampBeta = betaDf[defCols + [sampleName]].dropna().drop_duplicates().sort_values(['CpG_chrm', 'CpG_beg'])
    #sampBeta[sampleName] = [(beta + 10**-6) for beta in sampBeta[sampleName]]
    if mval == 'True':
        import math
        minNonZero = min([val for val in sampBeta[sampleName] if val > 0])
        maxNonOnes = max([val for val in sampBeta[sampleName] if val < 1])
        sampBeta.loc[sampBeta[sampleName] == 0, sampleName] = minNonZero
        sampBeta.loc[sampBeta[sampleName] == 1, sampleName] = maxNonOnes
        sampBeta[sampleName] = [math.log2(beta/(1 - beta)) for beta in sampBeta[sampleName]]
        return sampBeta
    else:
        return sampBeta

    
def convertBedToBigWig(bedFile, bwFileName, exePath, chrSizes):
    paras = [exePath, bedFile, chrSizes, bwFileName]
    os.system(' '.join(paras))
    return bwFileName    
    

def generateBedBigWig(betaFile, exePath, chrFile, outdir):
    """
    Generate bed and bigWig files from sorted beta value files
    """
    if outdir is None:
        basedir = os.path.basename(betaFile).split('.')[0]
        outdir = os.path.join(os.getcwd(), basedir)
        os.makedirs(outdir, exist_ok=True)
    else:
        outdir = os.path.join(outdir, basedir)
        os.makedirs(outdir, exist_ok=True)
    
    betaDf = pd.read_csv(betaFile, sep='\t')
    
    bedList = []
    bwList = []
    for sample in getSamples(betaDf):
        outBed = os.path.join(outdir, sample + '.bed')
        outBw = os.path.join(outdir, sample + '.bw')
        mDf = getMvalPerSample(betaDf, sample, mval='True')
        mDf.to_csv(outBed, sep='\t', index=False, header=False)
        convertBedToBigWig(outBed, outBw, exePath, chrFile)
        bedList.append(outBed)
        bwList.append(outBw)

        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate bed and bigWig files from sorted beta value files")
    parser.add_argument("-betaFile",
                    help="File with sorted beta values for the samples. ",
                    required=True)
    
    parser.add_argument("-exePath",
                    help="Path for bigWig exe file. ",
                    required=True)

    parser.add_argument("-chrFile",
                    help="Path for chr sizes file. ",
                    required=True)

    parser.add_argument("-outdir",
                    help="Output directory. Default: current directory. ",
                    required=False)
    args = parser.parse_args()
    betaFile = args.betaFile
    exePath = args.exePath
    chrFile = args.chrFile
    outdir = args.outdir
    generateBedBigWig(betaFile, exePath, chrFile, outdir)
