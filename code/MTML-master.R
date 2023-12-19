rm(list=ls())
library(dplyr)
library(MASS)
library(readxl)
library(mvtnorm)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(bigmemory)

source("./data/simu_data/simulate.R")
source("./code/single_trait.R")
source("./code/cauchyComb.R")
source("./code/ridge.R")
source("./code/MTML.R")
source("./code/statistic.R")

kk          <-   NULL                 
psmatrix    <-   NULL           
Likelihood  <-   "REML"       
svpal       <-   0.05              
Genformat   <-   1             
CLO         <-   16                   
snp         <-   4000   
spvalue     <-   0.05                  

#########################################################################################################
### Simulation

QTN              <-  read.csv("./data/simu_data/QTNs.csv",header = F)  # four QTNs
trait_info       <-  matrix(c(2, 5, 10, 500, 200, 100), 3, 2) 
genetic_effect1  <-  matrix(c(2.9161, 2.6082, 2.2588, 1.3041), nrow=4, ncol=1) 
genetic_effect2  <-  matrix(c(3.1944, 2.8571, 2.4743, 1.4286), nrow=4, ncol=1)

simu_data(QTN, trait_info, genetic_effect1, genetic_effect2)


##########################################################################################################
### MTML Method
         
 
genotype     <-  read.csv("./data/simu_data/snps.csv",header=T)
gen          <-  cbind(as.matrix(genotype[1:4000,1:3]),as.matrix(genotype[1:4000,4:202]))
trait_info   <-  matrix(c(2, 5, 10, 500, 200, 100), 3, 2)

MTML_perform <-  NULL
for (k in 1:nrow(trait_info)){
  for (r in 1:3){
    trait_k      <-  trait_info[k,1]
    N            <-  trait_info[k,2]
    fname        <-  paste("./data/simu_data/simu_phenotype/Y",r,"_",trait_k,".csv",sep="")
    phenotype    <-  read.csv(fname,header = T)
    phenotype    <-  as.matrix(phenotype[,-1])
    matrix_final <-  MTML(gen, phenotype, kk, psmatrix, N, svpal, Genformat, Likelihood, CLO, trait_k, snp)
    matrix_final[is.na(matrix_final)] = 9        
    statis       <-  statistic(gen, matrix_final, N, spvalue, snp, r, trait_k)
	MTML_perform <-  rbind(MTML_perform,statis)
  }
}
write.csv(MTML_perform,"./results/MTML_result.csv")
