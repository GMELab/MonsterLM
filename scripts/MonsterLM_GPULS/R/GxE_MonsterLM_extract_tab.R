#!/bin/bash
library(data.table)
library(MBESS)
args = (commandArgs(TRUE))
trait <- args[1]
raw_dir <- args[2]
results_dir <- args[3]
pheno_dir <- args[4]
setwd(raw_dir) #traits directory

gpuls_extract_tab_G_GxE = function(raw_dir) {
        
    setwd(raw_dir)

    GxE <- list.files(pattern = paste0(trait , ".*Gx.*"))
    G <- list.files(pattern = paste0(trait , ".*G\\..*"))
    
    x <- matrix(0,length(G)+length(GxE),1)
    df <- data.frame(trait=x,model=x,chr=x,split=x,pred_n=x,N=x,adj_r2=x,R2=x)

    for (i in 1:length(G)){
        #G results
        load(G[i])
        #print(G[i])
        pat1=sub("_G\\.RData.*","",G[i])
        pat2=sub(".*chr","",pat1)
        df[i,] <- c(trait,"G",sub("_.*","",pat2),sub(".*_","",pat2),length(contents$beta),length(contents$ypred),contents$adj_r2,contents$r2)
    }

    for (i in 1:length(GxE)){   
        #GxE results
        load(GxE[i])
        #print(GxE[i])
        pat1=sub("_G\\.RData.*","",G[i])
        pat2=sub(".*chr","",pat1)
        df[i+length(G),] <- c(trait,"GxE",sub("_.*","",pat2),sub(".*_","",pat2),length(contents$beta),length(contents$ypred),contents$adj_r2,contents$r2)
    }

    return(df)

}

sum_lmadj_G_GxE = function(df, pheno_dir) {

G_tot <- sum(as.numeric(df[df$model=="G",7]))
GxE_tot <- sum(as.numeric(df[df$model=="GxE",7]))

load(paste0(pheno_dir, "/phenotype/EonP_adjr2_" , trait , "_v1.1.RData"))
#load(paste0(pheno_dir, "/exposure/WHR_final_" , trait, "_v1.1.RData")); if we ever use >1 exposure this could be useful

x <- matrix(0,1,1)
df_sum <- data.frame(trait=x,G_adj_her=x,LCI_G_her=x,UCI_G_her=x,variance_G=x,GxE_adj_her=x,LCI_GxE_her=x,UCI_GxE_her=x,variance_GxE=x,E_her=x,LCI_E_her=x,UCI_E_her=x,variance_E=x,N=x,Pred_N_G=x,Pred_N_GxE=x,Pred_N_E=x,G_raw=x,GxE_raw=x)

df_sum[1,] <- c(trait, G_tot*(1-EonP), 0, 0, 0, GxE_tot*(1-EonP), 0, 0, 0, EonP, 0, 0, 0, df[1,6] , sum(as.numeric(df[df$model=="G",5])), sum(as.numeric(df[df$model=="GxE",5])), 1, G_tot, GxE_tot)

return(df_sum)

}


gpuls_extract_betas = function(raw_dir) {

    setwd(raw_dir)
    
    GxE <- list.files(pattern = paste0(trait , ".*Gx.*"))
    G <- list.files(pattern = paste0(trait , ".*G\\..*"))
    
    df <- matrix(0,0,6)
    colnames(df)[1:6] <- c("trait","model","chr","split","SNP_count","beta")
    m=function(x){as.matrix(x)}
    
    for (i in 1:length(G)){
        
        #G results
        load(G[i])
        #print(G[i])
        row_count=length(contents$beta)
        pat1=sub("_G\\.RData.*","",G[i])
        pat2=sub(".*chr","",pat1) #change chr to group for strat analysis
        df <- rbind(df, cbind(matrix(0,row_count,0) , m(rep(trait,row_count)),m(rep("G",row_count)),m(rep(sub("_.*","",pat2),row_count)),m(rep(sub(".*_","",pat2),row_count)),m(seq(1,row_count)),t(m(contents$beta))))

    }

    for (i in 1:length(GxE)){
        #GxE results
        load(GxE[i])
        #print(GxE[i])
        row_count=length(contents$beta)
        pat1=sub("_G\\.RData.*","",G[i])
        pat2=sub(".*chr","",pat1) #change chr to group for strat analysis
        df <- rbind(df, cbind(matrix(0,row_count,0) , m(rep(trait,row_count)),m(rep("GxE",row_count)),m(rep(sub("_.*","",pat2),row_count)),m(rep(sub(".*_","",pat2),row_count)),m(seq(1,row_count)),t(m(contents$beta))))

    }

    df_betas<-df
    return(df_betas)
}


mLM_adj_CI <- function(df, df_sum) {

vars<-c("variance","LCI","UCI")
df[,vars] <- c(NA,NA,NA)

    for (i in 1:nrow(df)) {

        n<-as.numeric(df[i,6])
        p<-as.numeric(df[i,5])
        adj_r2<-as.numeric(df[i,7])

    block_variance <- (1-as.numeric(df_sum$E_her)) * ( (n - 1)/(n - p - 1 ) )^2 * Variance.R2(adj_r2, n, p) 
    lci_block <- adj_r2 - 1.96*sqrt(block_variance)
    uci_block <- adj_r2 + 1.96*sqrt(block_variance)

    df[i,9:11] <- c(block_variance,lci_block,uci_block)

    }

    variance_total_G <- sum(as.numeric(df[df$model=="G",9]))
    variance_total_GxE <- sum(as.numeric(df[df$model=="GxE",9]))
    variance_total_E <- ( (as.numeric(df_sum$N) - 1)/(as.numeric(df_sum$N) - as.numeric(df_sum$Pred_N_E) - 1 ) ) * abs(Variance.R2(as.numeric(df_sum$E_her), as.numeric(df_sum$N), as.numeric(df_sum$Pred_N_E)))
    lci_total_G <- as.numeric(df_sum$G_adj_her) - 1.96*sqrt(variance_total_G)
    uci_total_G <- as.numeric(df_sum$G_adj_her) + 1.96*sqrt(variance_total_G)
    lci_total_GxE <- as.numeric(df_sum$GxE_adj_her) - 1.96*sqrt(variance_total_GxE)
    uci_total_GxE <- as.numeric(df_sum$GxE_adj_her) + 1.96*sqrt(variance_total_GxE)
    lci_total_E <- as.numeric(df_sum$E_her) - 1.96*sqrt(variance_total_E)
    uci_total_E <- as.numeric(df_sum$E_her) + 1.96*sqrt(variance_total_E)

df_sum[,c(3:5,7:9,11:13)] <- c(lci_total_G,uci_total_G,variance_total_G,lci_total_GxE,uci_total_GxE,variance_total_GxE,lci_total_E,uci_total_E,variance_total_E)

df_estimates_block_sum <- list("block_estimates"=df , "total_estimates"=df_sum)

return(df_estimates_block_sum)

} 

# MonsterLM Alogirthim to sum, adjust, and CI all blocks per trait
df <- gpuls_extract_tab_G_GxE(raw_dir) #all gpuls extracted parameters
df_sum <- sum_lmadj_G_GxE(df, pheno_dir) #all values summed
df_total <- mLM_adj_CI(df , df_sum) #all valus are adjsuted and CI calculated per block and per genome

write.csv(df_total$block_estimates, file = paste0(results_dir, trait, "_block_totals.csv"))
write.csv(df_total$total_estimates, file = paste0(results_dir, trait, "_total_estimate.csv"))

df_betas <- gpuls_extract_betas(raw_dir)
write.csv(df_betas, file = paste0(results_dir, trait, "_Pred_betas.csv"))