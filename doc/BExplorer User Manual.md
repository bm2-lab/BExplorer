# _BExplorer: optimizing base-editing gRNA designs and predicting pleiotropic effects_

## Introduction
BExplorer is an integrated and comprehensive computational pipeline for optimizing the design of gRNAs in silico

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
program requires 4 parameters, and each parameter is indispensable
    parser.add_argument('--chr', required=True, choices=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', 
        '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y'] , 
        help='input your target chromosome number(X, Y need to be capitalized)')
    parser.add_argument('--pos', required=True, type=int, help='input your target position number (GRCh38.p12)')
    parser.add_argument('--transtype', required=True, choices=['c2t', 'a2g'], help='c2t or a2t. \n '
                                                                                   '*****c2t(CBE): regard the base at your target position is C, '
                                                                                   'and you want to convert C to T. \n'
                                                                                   '******a2t(ABE):regard the base at your target position is A, and you want to convert A to G')
    parser.add_argument('--editor', required=True, choices= editor_list, help='choose one base editor.\n'
                                                                              '*****CBE: [ BE1, BE2, BE3, HF-BE3, BE4, BE4max, BE4-GAM, YE1-BE3, EE-BE3, YE2-BE3, YEE-BE3, '
                                                                              'VQR-BE3, VRER-BE3, SaBE3, SaBE4, SaBE4-Gam, SaKKH-BE3, Target-AID, BE-PLUS ]. \n'
                                                                              '*****ABE: [ ABE7.9, ABE7.10, xABE, ABESa, VQR-ABE, VRER-ABE, SaKKH-ABE ]')


* --chr
  Chromosome of target base. You can choose in 1~22, X, Y (Note: X, Y need to be capitalized)

* -g, --genome

  Reference genome file, must be hg19 or GRCh37 .


### example
1. git clone https://github.com/bm2-lab/BExplorer.git
2. cd BExplorer/src && tar -xzvf software.tar.gz 
3. install all required software and package
4. put `GRCh38.fna` into BExplorer/genomeDATA, put `dbsnp_146.hg38.vcf` into BExplorer/snpDATA
5. /BExplorer/python3 BExplorer.py --chr 1 --pos 155238215 --transtype c2t --editor BE3

For detailed usage information, please refer to the [BExplorer User Manual](/doc/BExplorer_User_Manual.md)

### Output

The output file "putative_neo.txt" contains all putative neoantigens information.

## Citation
Systematic exploration of optimized base-editing gRNA designs and pleiotropic effects with BExplorer, submitted, 2021

## Contact
zhanggc@tongji.edu.cn or qiliu@tongji.edu.cn

Tongji University, Shanghai, China

