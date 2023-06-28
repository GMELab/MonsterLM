#!/bin/bash
args = (commandArgs(TRUE))
trait <- c("ApoA" , "ApoB" , "Cholesterol" , "CRP" , "Glucose" , "HDL" , "HbA1c" , "height" , "LDL" , "TG" , "Total_bilirubin")#args[1]
raw_dir <- "/gmel/zerotide/common/GxE_MonsterLM/results/v1.1/"#args[2]
results_dir <-  "/gmel/zerotide/discipim/GxWHR_edits_R1/REAL_DATA/GxE_MonsterLM/results/v1.1/block_sums/"#args[3]
trait_length <- 11#(as.numeric(args[4]) - 1)
setwd(results_dir) #block_sums directory
#sum all mLM results
mLM_total<-as.numeric()
for(i in 1:trait_length){
mLM <- read.csv(file = paste0(results_dir , trait[i] , "_total_estimate.csv"))
mLM_total <- rbind(mLM_total , mLM[,-1])
}
write.csv(mLM_total , file = paste0(results_dir , "MonsterLM_Final_Results_v1.1.csv"))