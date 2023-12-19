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
svmlod <- 3                
Genformat <- 1             
CLO <- 8                   
N <- 3
trait_k <- 3

#****************************************************************************************************#
source("./compare_method_code/FASTmrEMMA/myFASTmrEMMA.R")

#****************************************************************************************************#
genotype <- read.csv("./data/real_data/gen_deal.csv",header=T)
gen <- cbind(as.matrix(genotype[,1:3]),as.matrix(genotype[,4:202]))
genRaw <- as.matrix(gen[,1:3])

phenotype <- read.csv("./data/real_data/trait_f3.csv",header = T)
phenotype<-as.matrix(phenotype[,4:ncol(phenotype)])
### add population structure
ph_vec_mat<-function(ph){
  ph<-as.matrix(ph)
  ph_matrix<-matrix(ph,199,trait_k)
  ph_matrix<-as.data.frame(ph_matrix)
  colnames(ph_matrix)<-paste("ph_",1:trait_k)
  return(ph_matrix=ph_matrix)
}
phen<-ph_vec_mat(phenotype)
phen <- as.matrix(phen)
str1<-read.table("./data/real_data/real_ph_f3.3.Q",header = F)
ps<-as.matrix(str1)

if (is.null(ps)==FALSE)
{
  ps1<-cbind(matrix(1,nrow=nrow(phen)),ps)
  vhat<-solve(crossprod(ps1,ps1))%*%crossprod(ps1,phen)
  vhat1<-vhat[-1,]
  phen_ph <- phen - ps%*%vhat1
}else{
  phen_ph<-phen 
}
phen_ph<-as.matrix(phen_ph) 

#****************************************************************************************************#
 
  phe <- as.matrix(phen_ph[,1])
  result <- myFASTmrEMMA(gen=gen,phe=phe,genRaw=genRaw,kk=NULL,psmatrix=NULL,svpal=0.05,svmlod=3,Genformat=1,Likelihood="REML",CLO=8)
  result1<-as.matrix(as.data.frame(result))
  
  phe <- as.matrix(phen_ph[,3])
  result <- myFASTmrEMMA(gen=gen,phe=phe,genRaw=genRaw,kk=NULL,psmatrix=NULL,svpal=0.05,svmlod=3,Genformat=1,Likelihood="REML",CLO=8)
  result3<-as.matrix(as.data.frame(result))
  
  result <- rbind(result1,result3)
  
  write.csv(result,"./results/FASTmrEMMA_real_result.csv")
