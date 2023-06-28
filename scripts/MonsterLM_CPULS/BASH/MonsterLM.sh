#!/bin/bash
Loc="/mnt/your/phenotypes/source_variables/" #traits directory
phen_out_dir="/mnt/your/results/Main_Analysis/step1_phenotypes_exposures/" #processed phenotype directory for Step1
out_dir="/mnt/your/results/Main_Analysis/step2_results/" #results directory for Step2 block fits
results_dir="/mnt/your/results/Main_Analysis/step3_results/" #results directory for Step3 estimate totals
script="/mnt/your/scripts/main/" # script directory

cd ${script}

#Gene-Environment Analysis (4 outcomes x 12 x 4 exposures for 192 outcome-exposure combinations)

i=1
for j in {1..4} # phenotype index
do
for cycle in {1..4} # exposure index
do
nohup Rscript ./MonsterLM.R $(($i + $((0 + 12 * $cycle)))) $j ${phen_out_dir} ${out_dir} ${results_dir} && PID1=$! &
nohup Rscript ./MonsterLM.R $(($i + $((1 + 12 * $cycle)))) $j ${phen_out_dir} ${out_dir} ${results_dir} && PID2=$! &
nohup Rscript ./MonsterLM.R $(($i + $((2 + 12 * $cycle)))) $j ${phen_out_dir} ${out_dir} ${results_dir} && PID3=$! &
nohup Rscript ./MonsterLM.R $(($i + $((3 + 12 * $cycle)))) $j ${phen_out_dir} ${out_dir} ${results_dir} && PID4=$! &
nohup Rscript ./MonsterLM.R $(($i + $((4 + 12 * $cycle)))) $j ${phen_out_dir} ${out_dir} ${results_dir} && PID5=$! &
nohup Rscript ./MonsterLM.R $(($i + $((5 + 12 * $cycle)))) $j ${phen_out_dir} ${out_dir} ${results_dir} && PID6=$! &
nohup Rscript ./MonsterLM.R $(($i + $((6 + 12 * $cycle)))) $j ${phen_out_dir} ${out_dir} ${results_dir} && PID7=$! &
nohup Rscript ./MonsterLM.R $(($i + $((7 + 12 * $cycle)))) $j ${phen_out_dir} ${out_dir} ${results_dir} && PID8=$! &
nohup Rscript ./MonsterLM.R $(($i + $((8 + 12 * $cycle)))) $j ${phen_out_dir} ${out_dir} ${results_dir} && PID9=$! &
nohup Rscript ./MonsterLM.R $(($i + $((9 + 12 * $cycle)))) $j ${phen_out_dir} ${out_dir} ${results_dir} && PID10=$! &
nohup Rscript ./MonsterLM.R $(($i + $((10 + 12 * $cycle)))) $j ${phen_out_dir} ${out_dir} ${results_dir} && PID11=$! &
nohup Rscript ./MonsterLM.R $(($i + $((11 + 12 * $cycle)))) $j ${phen_out_dir} ${out_dir} ${results_dir} && PID12=$! &
wait $PID1 $PID2 $PID3 $PID4 $PID5 $PID6 $PID7 $PID8 $PID9 $PID10 $PID11 $PID12
done # outcome loop
done # 

#######
