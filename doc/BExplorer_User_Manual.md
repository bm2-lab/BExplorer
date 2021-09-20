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
### Parameters

* --chr

  chromosome of target base. You can choose in `1~22, X, Y` (Note: X, Y need to be capitalized)

* --pos

  position number of target base(GRCh38.p12).
  
* --transtype

  select base conversion type. There are two options: `c2t and a2g`

* --editor

  base editor you want to use. You can choose 1 from 26 base editors.
  * CBE: `BE1, BE2, BE3, HF-BE3, BE4, BE4max, BE4-GAM, YE1-BE3, EE-BE3, YE2-BE3, YEE-BE3, VQR-BE3, VRER-BE3, SaBE3, SaBE4, SaBE4-Gam, SaKKH-BE3, Target-AID, BE-PLUS`
  * ABE: `ABE7.9, ABE7.10, xABE, ABESa, VQR-ABE, VRER-ABE, SaKKH-ABE`
  
> NOTE: if --transtype = c2t, --editor can only choose CBE. if --transtype = a2g, --editor can only choose ABE. 

### How to select parameters
1. Before selecting parameters, please confirm the chromosome number and position of your target site and the base editor you want to use.
2. Parameter --chr: please enter the target site chromosome after this parameter.
3. Parameter --pos: please enter the target site position after this parameter.
4. Parameter --transtype: please enter the base conversion type after this parameter. Please note that this parameter needs to be matched with the parameter --editor.
5. Parameter --editor: please enter base editor you want to use after this parameter. Please note that this parameter needs to be matched with the parameter --transtype.
6. Command line: /BExplorer/python3 BExplorer.py --chr \[target site chromosome] --pos \[target site position] --transtype \[base conversion type] --editor \[base editor you want to use]


### Example
1. git clone https://github.com/bm2-lab/BExplorer.git
2. cd BExplorer/src && tar -xzvf software.tar.gz 
3. install all required software and package
4. put `GRCh38.fna` into BExplorer/genomeDATA, put `dbsnp_146.hg38.vcf` into BExplorer/snpDATA
5. /BExplorer/python3 BExplorer.py --chr 1 --pos 155238215 --transtype c2t --editor BE3


### Output

The output file "output_file_`chr`:`pos`.tsv" contains all gRNA prediction information.

## Citation
Systematic exploration of optimized base-editing gRNA designs and pleiotropic effects with BExplorer, submitted, 2021

## Contact
zhanggc@tongji.edu.cn or qiliu@tongji.edu.cn

Tongji University, Shanghai, China

