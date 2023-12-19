rm(list=ls())
library(dplyr)
library(MASS)
library(readxl)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(bigmemory)
kk <- NULL                 
psmatrix <- NULL           
Likelihood <- "REML"       
svpal <- 0.05              
Genformat <- 1             
CLO <- 8                   
snp <- 215947   
spvalue<- 0.05   
           
source("./MTML_code_data/code/single_trait.R")
source("./MTML_code_data/code/cauchyComb.R")
source("./MTML_code_data/code/ridge.R")
source("./MTML_code_data/code/MTML.R")

genotype <- read.csv("./MTML_code_data/data/real_data/gen_deal.csv",header=T)
gen <- cbind(as.matrix(genotype[,1:3]),as.matrix(genotype[,4:202]))
trait_k <- 3
N <- 1
phenotype <- read.csv("./MTML_code_data/data/real_data/trait_f3.csv",header = T)
phenotype <- as.matrix(phenotype[,4:ncol(phenotype)])

### add population structure
ph_vec_mat <- function(ph){
  ph <- as.matrix(ph)
  ph_matrix <- matrix(ph,199,trait_k)
  ph_matrix <- as.data.frame(ph_matrix)
  colnames(ph_matrix) <- paste("ph_",1:trait_k)
  return(ph_matrix=ph_matrix)
}
phen <- ph_vec_mat(phenotype)
phen <- as.matrix(phen)
str1 <- read.table("./MTML_code_data/data/real_data/real_ph_f3.3.Q",header = F)
ps <- as.matrix(str1)

if (is.null(ps)==FALSE)
{
  ps1 <- cbind(matrix(1,nrow=nrow(phen)),ps)
  vhat <- solve(crossprod(ps1,ps1))%*%crossprod(ps1,phen)
  vhat1 <- vhat[-1,]
  phen_ph <- phen - ps%*%vhat1
}else{
  phen_ph <- phen 
}
phen_ph <- as.matrix(phen_ph) 

matrix_final <- MTML(gen,phen_ph,kk,psmatrix,N,svpal,Genformat,Likelihood,CLO,trait_k,snp)
matrix_final <- na.omit(matrix_final)
colnames(matrix_final) <- c("Id","Chromosome","Position","pvalue")
matrix_final <- data.frame(matrix_final)
matrix_final[,"pvalue"] <- as.numeric(matrix_final[,"pvalue"])
trait_final <- filter(matrix_final, pvalue < spvalue) 
write.csv(trait_final,"./MTML_code_data/results/MTML_real_result.csv")

 


