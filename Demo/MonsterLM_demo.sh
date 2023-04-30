#!/bin/bash
BASH_dir="~/MonsterLM/Demo" #set up the Demo directory from the .zip package

cd ${BASH_dir}

echo -e "
Scenario: Consider unadjusted values for 'outcome A' and 'exposure A' in 5,000 individuals.\n
Expected adjusted R^2 for Heritability: 0.2.
Expected adjusted R^2 for GxE: 0.05.\n
Run MonsterLM to perform the 3 steps outlined in the method overview to get estimates.
"

Rscript ./MonsterLM_demo.R "outcome_A.RData" "exposure_A.RData" "genotype_matrix_sim_3blocks.RData" ${BASH_dir}

echo -e "
Scenario: Consider unadjusted values for 'outcome B' and 'exposure B' in 5,000 individuals.\n
Expected adjusted R^2 for Heritability: 0.0.
Expected adjusted R^2 for GxE: 0.0.\n
Run MonsterLM to perform the 3 steps outlined in the method overview to get estimates.
"

Rscript ./MonsterLM_demo.R "outcome_B.RData" "exposure_B.RData" "genotype_matrix_sim_3blocks.RData" ${BASH_dir}
