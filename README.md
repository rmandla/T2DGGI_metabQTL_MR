# T2DGGI metabQTL Mendelian Randomization Scripts

This repository contains pipelines used to run MR for the T2DGGI GWAS summary statistics against mQTL datasets.

## Setup

`git clone https://github.com/rmandla/T2DGGI_metabQTL_MR.git`

### Requirements

* Python 3
  * pandas
  * numpy
* R (This code was tested on R 4.1)
  * TwoSampleMR
  * magrittr
  * tidyverse
  * vroom
  * MRPRESSO
* PLINK (This code was tested on PLINK1)

## Tutorial

1. [Genome-wide analyses](#gwa)
2. [Colocalizing SNPs only](#cso)

### Genome-wide analyses <a name="gwa"></a>

#### Data preparation

##### Run clumping of input GWAS summary statistics

```
python T2DGGI_metabQTL_MR/scripts/prep_MR.py clump \
  -s ../MR/EUR_MetalFixed_LDSC-CORR_Results1TBL-HG19-rsid-v2.gz \
  -r /humgen/florezlab/users/rmandla/1000G_hg19/1000G_hg19 \
  -o TEST1 \
  -e GWAS \
  -p plink
```

Here, clumping is run against input GWAS in HG19 using reference data from 1000G, where `-r` represents the path containing the reference files in PLINK bim/bed/fam format and `-s` represents the input summary statistics. `-o` represents the output header, such that the output files from the preparation will be `${output_header}.ALL.clumped`. `-e` represents the exposure type, and `-p` is the path to an executable PLINK binary. Since this step runs clumping on the input GWAS, and we are testing for MR between the GWAS and different metabolites, this step only needs to be run once for all genome-wide GWAS -> mQTL MR tests.

##### Run clumping of input mQTL summary statistics

```
python T2DGGI_metabQTL_MR/scripts/prep_MR.py clump \
  -s glucose_sst.txt.gz \
  -r /humgen/florezlab/users/rmandla/1000G_hg19/1000G_hg19 \
  -o TEST2 \
  -e mqtl \
  -p plink \
  -d UKB \
  -c CHROM \
  -b POS \
  -ea ALT \
  -nea REF \
  -pv logp \
  -isLogP True \
```

Here, clumping is run against an input mQTL dataset. The arguments are the same as above, with the `-s` now pointing to the path of the mQTL summary statistics. Additional arguments are also added to specify how to parse the mQTL dataset. Specifically, `-c` refers to the chromosome column name, `-b` refers to the base-pair column name, `-ea` refers to the effect allele column name, `-nea` refers to the non-effect allele column name, `-pv` refers to the p-value column name, and `-isLogP` specifies if the p-values in `-pv` are in -log10 format. `-d` refers to the input dataset, which is used to name the output file. Since this step clumps the mQTL data, it needs to be run for every metabolite tested. If there are no genome-wide significant hits in the mQTL data, no output file will be generated.

#### Running MR

##### Test for MR in the direction of GWAS -> mQTL

```
python T2DGGI_metabQTL_MR/scripts/prep_MR.py run_mr \
    -s ../MR/EUR_MetalFixed_LDSC-CORR_Results1TBL-HG19-rsid-v2.gz \
    -p glucose \
    -d UKB \
    -ex GWAS \
    -oc mqtl \
    -c CHROM \
    -b POS \
    -se se \
    -pv logp \
    -ea ALT \
    -nea REF \
    -isLogP True \
    -N 10000 \
    -a af \
    -be beta \
    -dpath glucose_sst.txt.gz \
    -cs TEST1.GWAS.T2DGGI.ALL.clumped \
    -o TEST1/TEST1
```

Here, we use the clumped output above to run MR only using the clumped SNPs. `-s` refers to the input GWAS summary statistic file. `-p` refers to the name of the metabolite, which is used in naming the output file. `-d` is the mQTL dataset name, used to name the output files. Additional flags also specify the column names for the required data from the mQTL, including `-c` for chromosome column, `-b` for position column, `-se` for standard error column, `-pv` for p-value column, `-ea` for effect allele column, `-nea` for non-effect allele column, `-isLogP` to specify if the p-values are in -log10 format, `-N` to specify either a column in the mQTL data with sample sizes OR an integer of the total sample size, `-af` for the allele frequency column, and `-be` for the beta column. `-ex` sets the exposure type and `-oc` sets the outcome type. These values must be GWAS or mqtl. `-dpath` is the path to the input mQTL summary statistic file. `-cs` is the path to the clumped SNPs generated in the clumping step above. And `-o` is the output header of the MR files.

##### Test for MR in the direction of pQTL -> GWAS

```
python T2DGGI_metabQTL_MR/scripts/prep_MR.py run_mr \
    -s ../MR/EUR_MetalFixed_LDSC-CORR_Results1TBL-HG19-rsid-v2.gz \
    -p glucose \
    -d ukb \
    -oc GWAS \
    -ex mqtl \
    -c CHROM \
    -b POS \
    -se se \
    -pv logp \
    -ea ALT \
    -nea REF \
    -isLogP True \
    -N 10000 \
    -a af \
    -be beta \
    -dpath glucose_sst.txt.gz \
    -cs TEST2.mQTL.UKB.ALL.clumped \
    -o TEST2/TEST2
```

The command for testing MR of mQTL -> GWAS is the same as for GWAS -> mQTL, except the exposure and outcomes are flipped, and the clumped SNPs reflect clumping on the mQTL summary statistics.

This script works by first reformatting and subsetting the input summary statistics to only include the clumped SNPs and all information necessary for MR. Then, it runs an Rscript to get the MR output.
