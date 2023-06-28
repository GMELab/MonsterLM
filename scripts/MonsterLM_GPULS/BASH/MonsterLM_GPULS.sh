#!/bin/bash
trait=("" "ApoA" "ApoB" "Cholesterol" "CRP" "Glucose" "HDL" "HbA1c" "height" "LDL" "TG" "Total_bilirubin") #sample traits list
tmp="/mnt/temp/data/" #temporary data directory to write to for GPULS; for treated phenotypes, exposures, genotypes
Loc="/mnt/your/source_data/Traits/" #traits directory
out_dir="/mnt/your/results_dir/" #results directory
out_dir2="/mnt/your/total/results/"
script="/mnt/zerotide/discipim/GxWHR_edits_R1/REAL_DATA/GxE_MonsterLM/scripts/"
cd ${script}

#PART 1#
#Run MonsterLM in RISTATP25 (Part 1; Treating Outcomes and Exposures)
for i in {1..11} #3 cycles of 4 traits; do in iterations of 4 traits
do
nohup Rscript ./GxE_MonsterLM_traits.R ${trait[i]} ${Loc} ${tmp} # treat phenotypes and exposures; will have to run this on ristatp25 for now
done #trait loop

#PART 2#
#Run MonsterLM in RISTATP servers (Part 2; Creating Trait Specific Genotypes)
i=1
for cycle in {0..2} #3 cycles of 4 traits; do in iterations of 4 traits
do
nohup Rscript ./GxE_MonsterLM_genotypes.R ${trait[$(($i + $((0 + 4 * $cycle))))]} ${Loc} ${tmp} && PID1=$! &
nohup Rscript ./GxE_MonsterLM_genotypes.R ${trait[$(($i + $((1 + 4 * $cycle))))]} ${Loc} ${tmp} && PID2=$! &
nohup Rscript ./GxE_MonsterLM_genotypes.R ${trait[$(($i + $((2 + 4 * $cycle))))]} ${Loc} ${tmp} && PID3=$! &
nohup Rscript ./GxE_MonsterLM_genotypes.R ${trait[$(($i + $((3 + 4 * $cycle))))]} ${Loc} ${tmp} && PID4=$! &
wait $PID1 $PID2 $PID3 $PID4
done #trait loop

#PART 3#
#Run in a server with GPULS available
cycle_START=1
cycle_END=4
for cycle in {0..2} #3 cycles of 4 traits; do in iterations of 4 traits
do
for group in {1..22} #chromosome loop
do
for block in {1..4} #split loop
do
for i in $(eval echo "{$((($cycle_START + $((4 * $cycle)))))..$((($cycle_END + $((4 * $cycle)))))}"); #trait loop
do
echo "Processing Phenotype ${trait[i]} on chromosome ${group} block ${block}; glhf"
# R2 for GxE only blocks
CUDA_VISIBLE_DEVICES=1,3,5,7,2,4,6,0 gpuls -x ${tmp}genotype/v1.1/UKB_09_${group}_${block}_1_${trait[i]}_G_STD.RData -y ${tmp}/phenotype/P_resid_${trait[i]}_v1.1.RData -a 350000 -b 36000 -d 10 -o ${out_dir2}/${trait[i]}_v1.1_chr${group}_${block}_G.RData && PID01=$! &
CUDA_VISIBLE_DEVICES=0,2,4,6,3,5,7,1 gpuls -x ${tmp}genotype/v1.1/UKB_09_${group}_${block}_1_${trait[i]}_GxWHR_QN.RData -y ${tmp}/phenotype/P_resid_${trait[i]}_v1.1.RData -a 350000 -b 36000 -d 10 -o ${out_dir2}/${trait[i]}_v1.1_chr${group}_${block}_GxWHR.RData && PID02=$! &
wait $PID01 $PID02
#sleep 16h && rm ${tmp}genotype/v1.1/UKB_09_${group}_${block}_1_${trait[i]}_* & #ideally place here but will not work because of permissions issues
done # trait loop
done # chromosome loop
done # block loop
done # cycle loop

#removal of blocks scripts; temporary until all permissions are sorted 
cycle_START=1
cycle_END=4
for cycle in {0..2} #3 cycles of 4 traits; do in iterations of 4 traits
do
for group in {1..22} #chromosome loop
do
echo "Removing all chr ${group} blocks from batch ${cycle} in 7 hours plus" && date
sleep 7h
for i in $(eval echo "{$((($cycle_START + $((4 * $cycle)))))..$((($cycle_END + $((4 * $cycle)))))}"); #trait loop
do
rm ${tmp}genotype/v1.1/UKB_09_${group}_*_1_${trait[i]}_*
echo "All chr ${group} blocks from batch ${cycle} have been removed as of" && date
done # trait loop
done # chromosome loop
done # cycle loop

#Run MonsterLM in RISTATP25 (Part 4; compiles all MonsterLM estimates)
for cycle in {0..2} #3 cycles of 4 traits; do in iterations of 4 traits
do
for i in $(eval echo "{$((($cycle_START + $((4 * $cycle)))))..$((($cycle_END + $((4 * $cycle)))))}"); #trait loop
do
nohup Rscript ./GxE_MonsterLM_extract_tab.R ${trait[i]} ${out_dir2} ${out_dir}/block_sums/ ${tmp} & #extract, sum, adjust, and CI results
done # trait loop
done # cycle loop
wait #
nohup Rscript ./GxE_MonsterLM_total_estimates.R ${trait} ${out_dir2} ${out_dir}/block_sums/ echo "${#trait[@]}" #compiles the final MonsterLM estimated for all traits
