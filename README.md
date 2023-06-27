![MonsterLM_logo 0 2](https://user-images.githubusercontent.com/80165657/235058821-108ef6a3-b631-4a67-9539-fcafc05dbd09.png)
# MonsterLM

A method to estimate the proportion of variance explained by heritability and gene-by-environment interactions, in a fast, accurate, efficient, and unbiased manner on biobank-scale datasets (N>300,000). 

# Table of Contents
- [Method Overview](#1)
- [Documentation](#paragraph1)
- [System Requirements and Installation Guide](#paragraph2)
- [Demo Requirements](#paragraph1)
- [Hardware Requirements](#paragraph2)
- [Software Requirements](#paragraph1)
- [Demo (with instructions for use)](#paragraph2)
- [License](#paragraph2)
  
## Method Overview <a name="1"></a>

The MonsterLM algorithm is a method designed to estimate the proportion of variance explained by additive genetic variants (h<sup>2</sup>) and gene-by-environment interactions (GxE) in a fast, accurate, efficient, and unbiased manner, particularly for large datasets with many single nucleotide polymorphisms (SNPs) compared to participants. The algorithm partitions autosomes into non-overlapping regions to estimate genome-wide heritability or interactions with environmental exposures, while minimizing residual population stratification effects and LD spills. The algorithm is based on a multiple linear regression modeling framework and can use CPUs with sufficient RAM support (>200 GB) or graphics processing units (GPUs) for computational speed gains.

See the figure below for an overview to apply MonsterLM:

![image](https://user-images.githubusercontent.com/80165657/235231670-a1e08a3b-7c2e-4da5-9459-6da963abdbd5.png)

## Documentation

MonsterLM is currently under peer review and official documentation will be posted upon acceptance.

## System Requirements and Installation Guide

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

## Demo (with instructions for use)

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

## License

GNU General Public License v3.0
