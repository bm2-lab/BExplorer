# _BExplorer: optimizing base-editing gRNA designs and predicting pleiotropic effects_

## Introduction
BExplorer is an integrated and comprehensive computational pipeline for optimizing the design of gRNAs in silico. BExplorer could design best gRNA for 26 existing base editors in base editing researching, and evaluate the pleiotropic effect of the corresponding base editing loci.

## Requirement
### software
* python3
* python2.7
* vcftools
* cas-offinder (need openCL)
* R
* CFD_Scoring (No user installation required)
> The software are packed in `src/software.tar.gz`
### python package
* pandas
* numpy
* regex
### R package
* RobustRankAggreg

## Usage
### input
1. Reference genome file: `GRCh38.fna`(fasta format, such as .fna, .fa, .fasta)
2. human reference SNPs file: `dbsnp_146.hg38.vcf` or other reference SNPs file base on hg38

### example
1. git clone https://github.com/bm2-lab/BExplorer.git
2. cd BExplorer/src && tar -xzvf software.tar.gz 
3. install all required software and package
4. put `GRCh38.fna` into BExplorer/genomeDATA, put `dbsnp_146.hg38.vcf` into BExplorer/snpDATA
5. /BExplorer/python3 BExplorer.py --chr 1 --pos 155238215 --transtype c2t --editor BE3

For detailed usage information, please refer to the [BExplorer User Manual](/doc/BExplorer_User_Manual.md)

## Citation
Systematic exploration of optimized base-editing gRNA designs and pleiotropic effects with BExplorer, submitted, 2021

## Contact
zhanggc@tongji.edu.cn or qiliu@tongji.edu.cn

Tongji University, Shanghai, China
