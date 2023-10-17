
R
source("plink.map.R")
source("plink.ped.R")
source("plink.ph.R")

genotype <- read.csv("gen_deal.csv",header=T)
genotype <- as.matrix(genotype[,4:202])
genotype[genotype==0] <- "11"
 genotype[genotype==1] <- "12"
 genotype[genotype==2] <- "22"
 gene <- genotype
 gene <- t(gene)
phenotype <- read.csv("trait_f3.csv",header = T)
phenotype<-as.matrix(phenotype[,4:ncol(phenotype)])
trait_k <- 3
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
str1<-read.table("real_ph_f3.3.Q",header = F)
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
phen <- phen_ph


n <- nrow(gene) 
q <- ncol(gene) 
p <- ncol(phen)
mapph <- plink.map(q)
write.table(mapph,"trait_f3.map",row.names=F,col.names=F,quote=F)
peddh <- plink.ped(n,phen,gene,p,q)
write.table(peddh,"trait_f3.ped",row.names=F,col.names=F,quote=F)
phh <- plink.ph(n,phen,p)
write.table(phh,"trait_f3.txt",row.names=F,col.names=F,quote=F)
q()



#****************************************** Plink bed et al ************************************************
./plink --file trait_f3 --make-bed --out trait_f3

#******************************************** Gemma related matrix************************************************
./gemma-0.98.5-linux-static-AMD64 -bfile trait_f3 -p trait_f3.txt -gk 1 -o relatedf3


# #******************************************* Gemma ************************************************
for ((i=1;i<=3;i++))
do
  ./gemma-0.98.5-linux-static-AMD64 -bfile trait_f3 -p trait_f3.txt  -k relatedf3.cXX.txt -lmm 1 -n $i -o gemmaf3$i
  
done

R
library(dplyr)
spvalue <- 0.00000023
gemma1 <- read.table("gemmaf31.assoc.txt",header=T)
gemma<- cbind(gemma1[,2],gemma1[,ncol(gemma1)])
colnames(gemma)<- c("Id","pvalue")
gemma <- data.frame(gemma)
gemma[,"pvalue"] <- as.numeric(gemma[,"pvalue"])
trait1 <- filter(gemma, pvalue < spvalue) 


gemma2 <- read.table("gemmaf32.assoc.txt",header=T)
gemma<- cbind(gemma2[,2],gemma2[,ncol(gemma1)])
colnames(gemma)<- c("Id","pvalue")
gemma <- data.frame(gemma)
gemma[,"pvalue"] <- as.numeric(gemma[,"pvalue"])
trait2 <- filter(gemma, pvalue < spvalue) 

gemma3 <- read.table("gemmaf33.assoc.txt",header=T)
gemma<- cbind(gemma3[,2],gemma3[,ncol(gemma1)])
colnames(gemma)<- c("Id","pvalue")
gemma <- data.frame(gemma)
gemma[,"pvalue"] <- as.numeric(gemma[,"pvalue"])
trait3 <- filter(gemma, pvalue < spvalue) 