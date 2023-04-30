#################################################################
## MonsterLM: Demonstration of Heritability and GxE Estimation ##
#################################################################

print("MonsterLM: Demonstration of Heritability and GxE Estimation")

#### All outcome, exposure, and genotype data is artifically 
#### simulated; no real-world values are used in this demo 

#Scenario: Consider raw values for 'outcome A' and 'exposure A' in 5,000 individuals
#Expected adjusted R^2 for Heritability: 0.2
#Expected adjusted R^2 for GxE: 0.05

#Scenario: Consider raw values for 'outcome B' and 'exposure B' in 5,000 individuals
#Expected adjusted R^2 for Heritability: 0.0
#Expected adjusted R^2 for GxE: 0.0

#### Step 0: Install BASH and R function dependencies for MonsterLM ####

#BASH arguments for simulated data
args = (commandArgs(TRUE))
outcome <- args[1]
exposure <- args[2]
genotype <- args[3]
setwd(args[4])

#R Library
library(MBESS)

#R Functions
loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}

quantNorm = function(x){qnorm(rank(x,ties.method = "average")/(length(x)+1))}

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

#################################################
### Step 1: Phenotype Outcome Quality Control ###
print("Step 1: Phenotype Outcome Quality Control")

outcome <- loadRData(outcome) # 5,000 x 1 matrix for raw values of dummy outcome A
exposure <- loadRData(exposure) # 5,000 x 1 matrix for raw values of dummy exposure A

P_final <-  quantNorm(outcome)
E_final <-  quantNorm(exposure)

P_resid <- quantNorm(resid(lm(P_final ~ E_final)))
EonP <- summary(lm(P_final ~ E_final))[[9]]

P_resid <- as.matrix(rm_Heteroscedasticity(P_resid, E_final))

######################################################################
##### Step 2: Run MonsterLM for an outcome-exposure combinations #####
print("Step 2: Run MonsterLM for an outcome-exposure combinations")

load(genotype) #load 5,000 x 1,500 block in additive ({0,1,2}) coding
genotype <- apply(genotype,2,standardization) #standardize genotype matrix

df_total <- as.data.frame(matrix(NA, nrow=1, ncol = 7))
colnames(df_total)[1:7] <- c("G | 0.2" , "GxE | 0.05" , "E | 0.0", "LCI_G", "UCI_G", "LCI_GxE", "UCI_GxE")

#Calculating Coefficients of Determination

for (o in 1:3) { #block loop for multiblock

  #split 5,000 X 1,500 simulated genotype block into 3 block of dimensions 5,000 x 500
  p = ncol(genotype)
  ranges = c(1,500,1000,1500)
    interaction_term = apply(genotype[,ranges[o]:ranges[o+1]] * E_final, 2, quantNorm)

	G_r2 = get_R2(as.matrix(genotype[,ranges[o]:ranges[o+1]]), P_resid)
	assign(paste0("G_r2_",o),G_r2[[1]][2])
	
	GxE_r2 = get_R2(interaction_term, P_resid)
	assign(paste0("GxE_r2_",o),GxE_r2[[1]][2])
	
}
   
##### Step 3: Sum, Adjust, and Confidence Intervals for Final Estimates  #####
print("Step 3: Sum, Adjust, and Confidence Intervals for Final Estimates")

df_total[,1] = sum(G_r2_1,G_r2_2,G_r2_3) * (1 - EonP) #sum 3 genotype blocks used for G estimation and adjust estimate by exposure variance
df_total[,2] = sum(GxE_r2_1,GxE_r2_2,GxE_r2_3) * (1 - EonP) #sum 3 genotype blocks used for GxE estimation and adjust estimate by exposure variance
df_total[,3] = EonP

#Confidence Intervals for G
var_tot=as.numeric()
for (k in 1:3) { #3 genotype blocks to calculate R^2 variance
adj_r2 = get(paste0("G_r2_", k))#
p = 500 # 500 SNP predictors per block (see Step 2 ranges split)
n = 5000 # 5,000 individuals per block
var_block = (1-EonP) * ( (n - 1)/(n - p - 1 ) )^2 * Variance.R2(adj_r2, n, p)
var_tot = rbind(var_tot , var_block)
}
df_total[,4:5] <- cbind(LCI = (df_total[,1] - 1.96*sqrt(sum(var_tot))), UCI = (df_total[,1] + 1.96*sqrt(sum(var_tot))))

#Confidence Intervals for GxE
var_tot=as.numeric()
for (k in 1:3) { #3 genotype blocks to calculate R^2 variance
adj_r2 = get(paste0("GxE_r2_", k))#
p = 500 # 500 SNP predictors per block (see Step 2 ranges split)
n = 5000 # 5,000 individuals per block
var_block = (1-EonP) * ( (n - 1)/(n - p - 1 ) )^2 * Variance.R2(adj_r2, n, p)
var_tot = rbind(var_tot , var_block)
}
df_total[,6:7] <- cbind(LCI = (df_total[,1] - 1.96*sqrt(sum(var_tot))), UCI = (df_total[,1] + 1.96*sqrt(sum(var_tot))))

print("Final Estimates with Confidence Intervals:")
df_total
