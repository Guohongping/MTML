MTML: An efficient multi-trait multi-locus GWAS method based on Cauchy combination test.
#' @param kk is the inital genomic relationship matrix.
#' @param psmatrix is the initial P-value.
#' @param Likelihood is the likelihood function type, "REML" represents restrictive maximum likelihood; "ML" represents maximum likelihood.
#' @param svpal is the P-value threshold for significance in the single trait part.
#' @param Genformat is the gene mode, including 0 and 1.
#' @param CLO is the number of computer cores.
#' @param trait_k is the insterest traits in GWAS, trait_k = d = (2,5,10).
#' @param snp is the markers being tested.
#' @param N represents execution/loop count, N = (500,200,100).
#' @param spvalue is the P-value threshold for significance in the multiple traits part.

# An simple example for MTML

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

source("simulate.R")
source("single_trait.R")
source("cauchyComb.R")
source("ridge.R")
source("MTML.R")
source("statistic.R")

kk <- NULL                 
psmatrix <- NULL           
Likelihood <- "REML"       
svpal <- 0.05              
Genformat <- 1             
CLO <- 8                   
snp <- 4000   
spvalue<- 0.05                  

#########################################################################################################
### Simulation

QTN <- read.csv("QTNs.csv",header = F)  # four QTNs
trait_info <- matrix(c(2,5,10,500,200,100),3,2) 
genetic_effect1 <- matrix(c(2.9161,2.6082,2.2588,1.3041),nrow=4,ncol=1) 
genetic_effect2 <- matrix(c(3.1944,2.8571,2.4743,1.4286),nrow=4,ncol=1)
simu_data(QTN,trait_info,genetic_effect1,genetic_effect2)


##########################################################################################################
### MTML Method
         
 
genotype <- read.csv("snps.csv",header=T)
gen <- cbind(as.matrix(genotype[1:4000,1:3]),as.matrix(genotype[1:4000,4:202]))
trait_info <- matrix(c(2,5,10,10,4,2),3,2)

MTML_perform <- NULL
for (k in 1:nrow(trait_info)){
  for (r in 1:3){
    trait_k <- trait_info[k,1]
    N <- trait_info[k,2]
    fname <- paste("Y",r,"_",trait_k,".csv",sep="")
    phenotype <- read.csv(fname,header = T)
    phenotype <- as.matrix(phenotype[,-1])
    matrix_final <- MTML(gen,phenotype,kk,psmatrix,N,svpal,Genformat,Likelihood,CLO,trait_k,snp)
    matrix_final[is.na(matrix_final)] = 9        
    statis <- statistic(gen,matrix_final,N,spvalue,snp,r,trait_k)
	MTML_perform <- rbind(MTML_perform,statis)
  }
}
