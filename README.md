![MonsterLM_logo 0 2](https://user-images.githubusercontent.com/80165657/235058821-108ef6a3-b631-4a67-9539-fcafc05dbd09.png)
# MonsterLM

A method to estimate the proportion of variance explained by heritability and gene-by-environment interactions, in a fast, accurate, efficient, and unbiased manner on biobank-scale datasets (N>300,000). 

# Table of Contents
- [Method Overview](#1)
- [System Requirements and Installation Guide](#3)
- [Demo (with instructions for use)](#4)
- [Running MonsterLM scripts on biobank-scale data](#5)
- [License](#6)
- [Contact Information](#7)
- [Citation](#8)
  
## Method Overview <a name="1"></a>

The MonsterLM algorithm is a method designed to estimate the proportion of variance explained by additive genetic variants (h<sup>2</sup>) and gene-by-environment interactions (GxE) in a fast, accurate, efficient, and unbiased manner, particularly for large datasets with many single nucleotide polymorphisms (SNPs) compared to participants. The algorithm partitions autosomes into non-overlapping regions to estimate genome-wide heritability or interactions with environmental exposures, while minimizing residual population stratification effects and LD spills. The algorithm is based on a multiple linear regression modeling framework and can use CPUs with sufficient RAM support (>200 GB) or graphics processing units (GPUs) for computational speed gains.

See the figure below for an overview to apply MonsterLM:

![image](https://user-images.githubusercontent.com/80165657/235231670-a1e08a3b-7c2e-4da5-9459-6da963abdbd5.png)

## System Requirements and Installation Guide <a name="3"></a>

MonsterLM can be run on all major platforms (e.g. GNU/Linux, macOS, Windows) with the R programming langage (R version 3.6.3 or newer).

### Demo Requirements

- Any standard computer (macOS, Linux, Windows).

### Hardware Requirements (Full-scale version)

- Requires a unix-like virtual environment supporting a minimum of 250GB RAM space for in-memory operations.
- (Optional) GPUs to run GPULS program for Step 2 Matrix inversion.

### Software Requirements

#### Essential Dependencies: programs
| Program | Description | Download |
| --- | --- | --- |
| BASH (≥ 5.0) | a unix shell and command language | [https://ubuntu.com/download/desktop](https://www.gnu.org/software/bash/) |
| R (≥ 4.0.0) or newer | R programming language | https://cran.r-project.org/ |
| GPULS (Optional) | Fast Ordinary Least Squares Computations using GPUs | https://gist.github.com/wjn0/fd1ded8a6e5033e5ca0d00ac131469ee |

#### Essential Dependencies: R packages
| R package | Install | Reference |
| --- | --- | --- |
| tidyverse | install.packages("tidyverse") | https://www.tidyverse.org/packages/ |
| data.table | install.packages("data.table") | https://cran.r-project.org/package=data.table |
| MBESS | install.packages("MBESS") | https://cran.r-project.org/web/packages/MBESS/index.html |
| gsl | install.packages("gsl") | https://cran.r-project.org/web/packages/gsl/index.html |

## Demo (with instructions for use) <a name="4"></a>

The following is a demonstration of the MosnterLM algorithim to estimate heritability (G) and gene-by-environment interactions (GxE). This analysis takes uses dummy data for 5,000 individuals and 1,500 SNPs. The full-scale analysis is designed to work on biobank scale data with >300,000 individuals and >1,000,000 SNPs. 

For this demo, all of the essential dependencies noted in the previous section are required except for GPULS, tidyverse, and data.table.

For demonstrative and efficiency purposes, the following analysis uses:

- a simulated genotype matrix of dimensions 5,000 individuals x 1,500 SNPs which gets partitioned into 3 matrices of 5,000 x 500.
- one non-null simulated outcome matrix (outcome_A) of dimensions 5,000 individuals x 1.
- one null simulated outcome matrix (outcome_B) of dimensions 5,000 individuals x 1.
- one non-null simulated exposure matrix (outcome_A) of dimensions 5,000 individuals x 1.
- one null simulated exposure matrix (outcome_B) of dimensions 5,000 individuals x 1.

This demo requires two scripts from the `MonsterLM/Demo` repository: i) the `MonsterLM_demo.sh` shell script and ii) the `MonsterLM_demo.R` R script. The shell script loops through the alogorithim estimations for non-null G and GxE (A) first and then for null G and GxE (B) estimations second. It uses steps 1, 2, and 3 in the R script that is visually depicted in the method overview. 

First, ensure that all dependencies are installed in your working environment. 

Next, set the BASH working directory in line 2 to your working directory
```
BASH_dir="~/MonsterLM/Demo" #set up the Demo directory from the .zip package
```

Once this is set up then grant the shell script executable permissions and run as follows

```
chmod +x MosnterLM_demo.sh
./MosnterLM_demo.sh
```

A succesful output will compute for the combinations outcome_A-exposure_A and outcome_B-exposure_B by running the full algorithim. The run time for this script on a standard computer should be <10 seconds. The output window should say:

```

Scenario: Consider unadjusted values for 'outcome A' and 'exposure A' in 5,000 individuals.

Expected adjusted R^2 for Heritability: 0.2.
Expected adjusted R^2 for GxE: 0.05.

Run MonsterLM to perform the 3 steps outlined in the method overview to get estimates.

[1] "MonsterLM: Demonstration of Heritability and GxE Estimation"
[1] "Step 1: Phenotype Outcome Quality Control"
[1] "Step 2: Run MonsterLM for an outcome-exposure combinations"
[1] "Step 3: Sum, Adjust, and Confidence Intervals for Final Estimates"
[1] "Final Estimates with Confidence Intervals:"
    G | 0.2 GxE | 0.05     E | 0.0     LCI_G     UCI_G    LCI_GxE    UCI_GxE
1 0.2025808 0.05126036 0.001621415 0.1711931 0.2339685 0.02586605 0.07665466

Scenario: Consider unadjusted values for 'outcome B' and 'exposure B' in 5,000 individuals.

Expected adjusted R^2 for Heritability: 0.0.
Expected adjusted R^2 for GxE: 0.0.

Run MonsterLM to perform the 3 steps outlined in the method overview to get estimates.

[1] "MonsterLM: Demonstration of Heritability and GxE Estimation"
[1] "Step 1: Phenotype Outcome Quality Control"
[1] "Step 2: Run MonsterLM for an outcome-exposure combinations"
[1] "Step 3: Sum, Adjust, and Confidence Intervals for Final Estimates"
[1] "Final Estimates with Confidence Intervals:"
      G | 0.2  GxE | 0.05       E | 0.0      LCI_G      UCI_G     LCI_GxE
1 0.007945522 0.005875249 -0.0001499552 -0.0151448 0.03103584 -0.01707937
     UCI_GxE
1 0.02882987


```

Specific coding techniques used for each step can be viewed for Steps 1 - 3 in `MonsterLM_demo.R`.

## Running MonsterLM on biobank-scale data <a name="5"></a>

Prior to running MonsterLM on biobank-scale data, first ensure individual-level genotype and phenotype data is prepared in the .bed format to be accessed by PLINK. Individal and SNP quality control is described in 4 steps in `MonsterLM/scripts/PLINK`:

  0. **`PLINK_QC.sh`**

Source .bed files should be pre-processed to meet the following criteria mentioned in Step 0.

```
#--------------------------------------------------------------------------------
#                Prepare genotypes for h2 and GxE testing
#--------------------------------------------------------------------------------
# Step 0: Preliminary QC for biobank data in .fam, .bim, and .bed including:
#
# Individual Exclusion criteria: (1) non-white related British ancestry, (2) high 
# ancestry-specific heterozygosity, (3) high genotype missingness (>0.05), 
# (3) mismatching genetic ancestry, (4) sex chromosome aneuploidy, 
# (5) mismatching gender sex and genetic sex, and (6) consent withdrawal 
# at the time of analysis.
#
# SNP exclusion criteria included: (1) SNPs with low imputation quality 
# (INFO score < 0.30), (2) call rate < 0.95, and (3) ambiguous or duplicated SNPs.
#
```

PLINK scripts can be set up according to user preference. The following offers and example template to prepare the genotypic and phenotypic files for MonsterLM.

Step 1: Define inputs and paths.

```
# Define inputs

maf=$1

# Define paths (dummy paths used below)

outdir="/your/h2/output_directory"
subDir="/your/reference/allele/directory" #
UKB_genotypes="/UKBIOBANK/username/UKB_genotypes_plink" #.bed source files location
plink19="/home/programs/plink_1.9/plink"
chr=$(echo 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
```

Step 2: SNP Quality Control.

```
# Script start


if [ $# = 1 ]; then


  #-------------------------------------------------------------------------------
  # Step 2: SNP quality control
  #-------------------------------------------------------------------------------
  # Hardy-Weingberg eqilibrium: 1e-10
  # MAF: 0.05
  # genotype missingness: 0.05
  #-------------------------------------------------------------------------------


  function hq_SNPs {

    set -e

    mkdir -p ${outdir}/plink19_UKB_MAF_${maf}

    for i in ${chr[@]}; do

      ${plink19} \
        \
        --noweb \
        \
        --bfile ${UKB_genotypes}/chr.${i}.UKB_v3_INFO_0.30.UNRELATED.BRITISH.RELATED.MAF_0.001.HWE.FINAL \
        \
        --maf ${maf} --geno 0.05 --hwe 1e-10 \
        \
        --write-snplist \
        \
        --out ${outdir}/plink19_UKB_MAF_${maf}/chr.${i} &

    done
    wait

    for i in ${chr[@]}; do

      ${plink19} \
        \
        --noweb \
        \
        --bfile ${UKB_genotypes}/chr.${i}.UKB_v3_INFO_0.30.UNRELATED.BRITISH.RELATED.MAF_0.001.HWE.FINAL \
        \
        --keep-allele-order \
        \
        --extract ${outdir}/plink19_UKB_MAF_${maf}/chr.${i}.snplist \
        \
        --make-bed \
        \
        --out ${outdir}/plink19_UKB_MAF_${maf}/chr.${i}.UKB_v3_INFO_0.30.UNRELATED.BRITISH.RELATED.MAF_${maf}.GENO_0.05.HWE.FINAL &

    done
    wait

  }
  hq_SNPs
```

Step 3: LD-prune Genotype Matrices and Recode to the Additive Model Format

```
  #-----------------------------------------------------------------------------
  # Step 3: LD-prune Genotype Matrices and Recode to Additive Model ({0,1,2})
  #-----------------------------------------------------------------------------


  function LD_prune_for_MonsterLM {

    set -e

    for i in ${chr[@]}; do

      ${plink19} \
        \
        --noweb \
        \
        --bfile ${outdir}/plink19_UKB_MAF_${maf}/chr.${i}.UKB_v3_INFO_0.30.UNRELATED.BRITISH.RELATED.MAF_${maf}.GENO_0.05.HWE.FINAL \
        \
        --keep-allele-order \
        \
        --indep-pairwise 1000 500 0.9 \
        \
        --out ${outdir}/plink19_UKB_MAF_${maf}/chr.${i} &

    done
    wait

    for i in ${chr[@]}; do

      ${plink19} \
        \
        --noweb \
        \
        --bfile ${outdir}/plink19_UKB_MAF_${maf}/chr.${i}.UKB_v3_INFO_0.30.UNRELATED.BRITISH.RELATED.MAF_${maf}.GENO_0.05.HWE.FINAL \
        \
        --keep-allele-order \
        \
        --extract ${outdir}/plink19_UKB_MAF_${maf}/chr.${i}.prune.in \
        \
        --make-bed \
        \
        --out ${outdir}/plink19_UKB_MAF_${maf}/chr.${i}.UKB_v3_INFO_0.30.UNRELATED.BRITISH.RELATED.MAF_${maf}.GENO_0.05.HWE.LD_0.9_for_MonsterLM &

    done
    wait

  }
  LD_prune_for_MonsterLM
  
function recode_for_MonsterLM {

    set -e

    for i in ${chr[@]}; do

      ${plink19} \
        \
        --noweb \
        \
        --bfile ${outdir}/plink19_UKB_MAF_${maf}/chr.${i}.UKB_v3_INFO_0.30.UNRELATED.BRITISH.RELATED.MAF_${maf}.GENO_0.05.HWE.LD_0.9_for_MonsterLM \
        \
        --recodeA \
        \
        --recode-allele ${subDir}/ref_allele_09_${chr}
        \
        --out ${outdir}/plink19_UKB_MAF_${maf}/chr.${i}.UKB_v3_INFO_0.30.UNRELATED.BRITISH.RELATED.MAF_${maf}.GENO_0.05.HWE.LD_0.9_for_MonsterLM_ADDITIVE &

    done
    wait
  }
  recode_for_MonsterLM

```

Step 4: Clear working directories of temporary PLINK files.

```
#-------------------------------------------------------------------------------
  # Step 4: Clear out
  #-------------------------------------------------------------------------------
  # Remove all temporary files
  #-------------------------------------------------------------------------------


  function clear_out {

    rm -rf ${outdir}/plink19_UKB_MAF_${maf}/*log
    rm -rf ${outdir}/plink19_UKB_MAF_${maf}/*frq
    rm -rf ${outdir}/plink19_UKB_MAF_${maf}/*prune*
    rm -rf ${outdir}/plink19_UKB_MAF_${maf}/*snplist

  }
  clear_out

else

  echo "ERROR: incorrect number of input arguments"

fi
```


## License <a name="6"></a>

GNU General Public License v3.0

## Contact Information <a name="7"></a>

Any queries pertaining to the MonsterLM scripts or methodological framework can be addressed to either: Matteo Di Scipio (discipim@mcmaster.ca) or Guillaume Paré (pareg@mcmaster.ca).

## Citation <a name="8"></a>

### *Nature Communications*

Provisionally accepted.

### Github repository

[ZENODO LINK]
