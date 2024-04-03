#!/bin/sh
library(data.table)
library(MBESS)
options(scipen=999)
args = (commandArgs(TRUE))
#
trait <- as.numeric(args[1])#
exposure_var <- as.numeric(args[2])#
Loc <- "/mnt/your/source_exposures_phenotypes/" # exposures and phenotypes directory
Loc2 <- "/mnt/your/source_covariates/" # covariates directory
phen_out_dir <- args[3]
out_dir <- args[4] #blocks output directory; MonsterLM step 2 results
results_dir <- args[5] #final results; MonsterLM step 3 results
env_type <- args[6] #"continuous" or "dichotomous"

# block specs
chromosome <- seq(1:22) #eval(parse(text=args[2]))
blocks_dir="/mnt/your/standardized_genotype_matrices/"
max_sets = c(4,4,4,3,3,4,4,4,3,3,3,3,2,2,2,2,2,2,2,2,1,1)

quantNorm = function(x)
{
	qnorm(rank(x,ties.method = "average")/(length(x)+1))
}

standardization = function(x)
{
	return((x - mean(x))/sd(x))
}

get_R2 = function(Datas, outcome)
{
	Data = list()
	n = dim(Datas)[1]
	m = dim(Datas)[2]

	C_all = t(Datas) %*% outcome	
	
	C_matrix = crossprod((Datas))
	inv_matrix = solve(C_matrix)
	rm(C_matrix)
	
	Betas = inv_matrix %*% C_all

	Predicted = Datas %*% Betas
	
	R2 = var(Predicted)
	adj_R2 = 1 - ( 1 - R2 ) * ( n - 1 ) / ( n - m - 1 )

	rm(C_all)
	rm(Datas)
	rm(inv_matrix)
	gc()
	
	Data[[1]] = c(R2, adj_R2)
	Data[[2]] = Betas
	
	return(Data)
}

rm_stratificaiton = function(EP)
{
	EP_final <- quantNorm(EP)
	EP_final <- resid(lm( EP_final ~ Age[, 3] + Sex[, 3] + PCs[, 3:22] ))
	EP_final <- quantNorm(EP_final)
	
	return(EP_final)
}

rm_Heteroscedasticity = function(P_resid, E_final)
{
	nb.partitions <- 20

	partition.size <- 1 / nb.partitions / 2
	E_ranks <- rank(E_final,ties.method = "average")/(length(E_final)+1)
	P_resid_final <- P_resid
	for( i in 1:nb.partitions) {
		lower_lower_bound <- ( i - 1 ) * partition.size
		lower_upper_bound <- i * partition.size
		upper_lower_bound <- 1 - i * partition.size
		upper_upper_bound <- 1 - ( i - 1 ) * partition.size
		select <- which( (E_ranks >=  lower_lower_bound & E_ranks <  lower_upper_bound ) | ( E_ranks > upper_lower_bound & E_ranks <= upper_upper_bound ) )
		P_resid_final[select] <-  standardization( P_resid_final[select] )
	}
	
	return(P_resid_final)
}

write_beta = function(bim, beta, id, tag)
{
	Data = cbind(bim[,2], beta)
	colnames(Data) = c("SNP", "beta")
	write.table(Data, file=paste0("betas/betas_chr_", chr, "_", id, "_", tag, ".txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
	
	return(0)
}

loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}

# Main function
setwd(Loc)

#load data
#set.seed(1)
print(paste0("trait number ",trait))
load( paste0(Loc , "outcomes_N.RData") ) # outcomes matrix 
trait = colnames(outcomes)[-1][trait]
print(paste0("trait name ",trait))
pheno = as.matrix( outcomes[ , which( colnames(outcomes) == trait ) ] )
rm(outcomes)

print(paste0("exposure number ",exposure_var))
load( paste0(Loc , "exposure_exposures_325989.RData") )
exposure_var = colnames(exposure_exposures)[-1][exposure_var]
print(paste0("exposure name ",exposure_var))
exposure = as.matrix( exposure_exposures[ , which( colnames(exposure_exposures) == exposure_var ) ] )
rm(exposure_exposures)
Age = as.matrix(fread(paste0(Loc2 , "/Age_orig.txt"), header = T))
Sex = as.matrix(fread(paste0(Loc2 , "/Sex_orig.txt"), header = T))
PCs = as.matrix(fread(paste0(Loc2 , "/PCs_orig.txt"), header = T))

#identify NAs (there is no NA in exposure, Age, Sex and PCs)
index = unique(c(which(is.na(pheno)), which(is.na(exposure))))
pheno = as.matrix(pheno)
exposure = as.matrix(exposure)
if(length(index) > 0)
{
	pheno = as.matrix(pheno[-index, ])
	exposure = as.matrix(exposure[-index, ])
	Age = Age[-index, ]
	Sex = Sex[-index, ]
	PCs = PCs[-index, ]
}

if (env_type == "continuous") {
E_final = rm_stratificaiton(exposure[,1]) #if E is dichotomous, skip
} else {
E_final = standardization(exposure[,1])
}
	
P_final = rm_stratificaiton(pheno[,1])
EonP <- summary(lm(P_final ~ E_final))[[9]]

P_resid <- resid(lm(P_final~ E_final))
P_resid <- quantNorm(P_resid)
if (env_type == "continuous") {
P_resid <- rm_Heteroscedasticity(P_resid, E_final) #if E is dichotomous, skip
} else {
    unique_Es = unique(E_final)
     for (unique_E_val in unique_Es) {
          P_resid[E_final == unique_E_val] <- quantNorm(P_resid[E_final == unique_E_val])
        }
}
P_resid <- as.matrix(P_resid)

save(EonP, file=paste0(phen_out_dir,"EonP_", exposure_var , "_on_", trait ,".RData"))
save(P_resid, file=paste0(phen_out_dir,"Presid_", exposure_var , "_in_", trait ,".RData"))
print(paste0("MonsterLM Gene-exposure Step 1: Outcome ", trait, " QC Complete; Exposure ", exposure_var, " Complete"))

top_blocks_G_r2 = c("chr", "set", "nb_indi", "nb_SNPs", "R2", "adj_R2")
top_blocks_GxE_only_r2 = c("chr", "set", "nb_indi", "nb_SNPs", "R2", "adj_R2")

setwd(out_dir) 
for (chr in 1:length(chromosome)) #chromosome we are on
{
for (set in 1:max_sets[chr]) #set per chromosome we are on
{
	load(paste0(blocks_dir, "/UKB_stdnorm_09_" , chr , "_" , set , "_1.RData"))
	if(length(index) > 0) #Remove the NA participants from genotype data
	{
		geno_data = apply(geno_data[-index, ], 2, standardization)
	}

	if (env_type == "continuous") {
	interaction_term = apply(geno_data*E_final, 2, quantNorm) # for continuous E
        } else {
        interaction_term = apply(geno_data*E_final, 2, standardization) # for dichotomous E
	}
		
	#GxE only
	GxE_only_r2 = get_R2(interaction_term, P_resid)
	item_GxE_only = c(chr, set, dim(interaction_term)[1], dim(interaction_term)[2], GxE_only_r2[[1]])			
	top_blocks_GxE_only_r2 = rbind(top_blocks_GxE_only_r2, item_GxE_only)
		
	rm(interaction_term)
	gc()
		
	#G only
	G_r2 = get_R2(geno_data, P_resid)
	item_G = c(chr, set, dim(geno_data)[1], dim(geno_data)[2], G_r2[[1]])			
	top_blocks_G_r2 = rbind(top_blocks_G_r2, item_G)
	
    print(paste0("MonsterLM Gene-exposure Step 2: Outcome ", trait, " and Exposure ", exposure_var, " Matrix Inversion Complete for chr ", chr, " set ", set, " is done"))

	rm(geno_data)
	gc()
}

setwd(out_dir)
write.table(top_blocks_G_r2, file=paste0("Gene-exposure_",trait," with ", exposure_var ,"_chr_", chr, "_G.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t")
write.table(top_blocks_GxE_only_r2, file=paste0("Gene-exposure_",trait," with ", exposure_var ,"_chr_", chr, "_GxE_only.txt"), col.names=F, row.names=F, quote=F, sep="\t")

}

headers <- c(top_blocks_G_r2[1,])
top_blocks_G_r2 <- as.data.frame(top_blocks_G_r2[-1,])
colnames(top_blocks_G_r2)[1:length(headers)] <- headers
top_blocks_G_r2 <- as.data.frame(apply(top_blocks_G_r2,2,as.numeric))

headers <- c(top_blocks_GxE_only_r2[1,])
top_blocks_GxE_only_r2 <- as.data.frame(top_blocks_GxE_only_r2[-1,])
colnames(top_blocks_GxE_only_r2)[1:length(headers)] <- headers
top_blocks_GxE_only_r2 <- as.data.frame(apply(top_blocks_GxE_only_r2,2,as.numeric))

#MonsterLM Step 3: Block Summation, Adjustments, and Confidence Intervals
setwd(results_dir)
#G total
df_total_G <- data.frame(model="G",outcome=trait,exposure=exposure_var,pred_N=sum(top_blocks_G_r2$nb_SNPs),N=top_blocks_G_r2$nb_indi[1],est_adj=(1-EonP)*sum(top_blocks_G_r2$adj_R2))

#Confidence Intervals for G
var_tot=as.numeric()
for (k in 1:nrow(top_blocks_G_r2)) {
adj_r2 = top_blocks_G_r2$adj_R2[k]
p = top_blocks_G_r2$nb_SNPs[k]
n = top_blocks_G_r2$nb_indi[k]
var_block = (1-EonP) * ( (n - 1)/(n - p - 1 ) )^2 * Variance.R2(adj_r2, n, p)
var_tot = rbind(var_tot , var_block)
}
df_total_G <- cbind(df_total_G,LCI = (df_total_G$est_adj - 1.96*sqrt(sum(var_tot))),UCI = (df_total_G$est_adj + 1.96*sqrt(sum(var_tot))),variance_total=sum(var_tot),standard_deviation=sqrt(sum(var_tot)))


#GxE_total
df_total_GxE_only <- data.frame(model="GxE",outcome=trait,exposure=exposure_var,pred_N=sum(top_blocks_GxE_only_r2$nb_SNPs),N=top_blocks_G_r2$nb_indi[1],est_adj=(1-EonP)*sum(top_blocks_GxE_only_r2$adj_R2))

#Confidence Intervals for GxE
var_tot=as.numeric()
for (k in 1:nrow(top_blocks_GxE_only_r2)) {
adj_r2 = top_blocks_GxE_only_r2$adj_R2[k]
p = top_blocks_GxE_only_r2$nb_SNPs[k]
n = top_blocks_GxE_only_r2$nb_indi[k]
var_block = (1-EonP) * ( (n - 1)/(n - p - 1 ) )^2 * Variance.R2(adj_r2, n, p)
var_tot = rbind(var_tot , var_block)
}
df_total_GxE_only <- cbind(df_total_GxE_only,LCI = (df_total_GxE_only$est_adj - 1.96*sqrt(sum(var_tot))),UCI = (df_total_GxE_only$est_adj + 1.96*sqrt(sum(var_tot))),variance_total=sum(var_tot),standard_deviation=sqrt(sum(var_tot)))

df_total_final=rbind(df_total_G,df_total_GxE_only)

print(paste0("MonsterLM Gene-exposure Step 3: Outcome ", trait, " and Exposure ", exposure_var, " Summation, Adjusments, and Confidence Intervals Formation Complete"))

write.table(df_total_final , file=paste0("Gene-exposure_",trait," with ", exposure_var ,"_total.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
