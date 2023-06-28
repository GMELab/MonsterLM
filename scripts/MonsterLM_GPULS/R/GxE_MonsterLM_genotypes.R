#!/bin/bash
library("data.table")
blocks_dir="/gmel/your/standardized/additive/blocks/UKB_v3/"
max_sets = c(4,4,4,3,3,4,4,4,3,3,3,3,2,2,2,2,2,2,2,2,1,1)
standardization = function(x)
{
	return((x - mean(x))/sd(x))
}
quantNorm = function(x)
{
	qnorm(rank(x,ties.method = "average")/(length(x)+1))
}
args = (commandArgs(TRUE))
trait <- args[1]
Loc <- args[2]
tmp <- args[3]
setwd(tmp) #traits directory
#load data
pheno = as.matrix(fread(paste0(Loc , trait, "/Pheno/" , trait, "_orig.txt"), header = T))
#identify NAs (there is no NA in WHR, Age, Sex and PCs)
index = which(is.na(pheno[,3]))

#G only blocks and GxE only blocks
for (i in 1:22){
  for (j in 1:max_sets[i]) {
  load(paste0(blocks_dir, "/UKB_stdnorm_09_" , i , "_" , j , "_1.RData"))
  if (length(index) > 0)
  {
  geno_data <- apply(geno_data[-index, ], 2, standardization)
  }
  save(geno_data, file=paste0(tmp, "genotype/v1.1/UKB_09_" , i , "_" , j , "_1_" , trait , "_G_STD.RData"))
  load(paste0(tmp, "/exposure/WHR_final_" , trait, "_v1.1.RData"))
  GxE <- apply(geno_data*as.vector(E_final),2,quantNorm)
  rm(geno_data)
  save(GxE, file=paste0(tmp, "genotype/v1.1/UKB_09_" , i , "_" , j , "_1_" , trait , "_GxWHR_QN.RData"))
  rm(GxE)
  }
}
