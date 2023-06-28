#!/bin/sh
library("data.table")
options(scipen=999)
#
args = (commandArgs(TRUE))
trait <- args[1]
Loc <- args[2]
tmp <- args[3]

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
	
	return(c(R2, adj_R2))
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


# Main function
setwd(Loc) #traits directory

block_ids = c(4, 4, 4, 3, 3, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1)

#load data
pheno = as.matrix(fread(paste0(Loc , trait, "/Pheno/" , trait, "_orig.txt"), header = T))
WHR = as.matrix(fread(paste0(Loc , "/WHR/Pheno/WHR_orig.txt"), header = T))
Age = as.matrix(fread(paste0(Loc , "/Age_orig.txt"), header = T))
Sex = as.matrix(fread(paste0(Loc , "/Sex_orig.txt"), header = T))
PCs = as.matrix(fread(paste0(Loc , "/PCs_orig.txt"), header = T))

#identify NAs (there is no NA in WHR, Age, Sex and PCs)
index = which(is.na(pheno[,3]))
if(length(index) > 0)
{
	pheno = pheno[-index, ]
	WHR = WHR[-index, ]
	Age = Age[-index, ]
	Sex = Sex[-index, ]
	PCs = PCs[-index,]
}

E_final = rm_stratificaiton(WHR[, 3])
P_final = rm_stratificaiton(pheno[, 3])

P_resid <- resid(lm(P_final~ E_final))
EonP <- summary(lm(P_final~ E_final))[[9]]
P_resid <- quantNorm(P_resid)

P_resid <- rm_Heteroscedasticity(P_resid, E_final)
P_resid <- as.matrix(P_resid)

save(EonP, file=paste0(tmp, "/phenotype/EonP_adjr2_" , trait, "_v1.1.RData"))
save(P_resid, file=paste0(tmp, "/phenotype/P_resid_" , trait, "_v1.1.RData"))
save(E_final, file=paste0(tmp, "/exposure/WHR_final_" , trait, "_v1.1.RData"))
