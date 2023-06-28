#!/bin/bash

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
#--------------------------------------------------------------------------------
# Step 1: Define Inputs and Paths
#-------------------------------------------------------------------------------- 

# Define inputs

maf=$1

# Define paths (dummy paths used below)

outdir="/your/h2/output_directory"
subDir="/your/reference/allele/directory" #
UKB_genotypes="/UKBIOBANK/username/UKB_genotypes_plink" #.bed source files location
plink19="/home/programs/plink_1.9/plink"
chr=$(echo 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)


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
