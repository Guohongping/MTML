

  ./gemma-0.98.5-linux-static-AMD64 -bfile trait_f3 -p trait_f3.txt  -k relatedf3.cXX.txt  -lmm 1 -n 1 2 3 -o mvlmmf3 

R
library(dplyr)
spvalue <- 0.1
mvlmm <- read.table("mvlmmf3.assoc.txt",header=T) 
mvlmm<- cbind(mvlmm[,2],mvlmm[,ncol(mvlmm)])
colnames(mvlmm)<- c("Id","pvalue")
mvlmm <- data.frame(mvlmm)
mvlmm[,"pvalue"] <- as.numeric(mvlmm[,"pvalue"])
trait <- filter(mvlmm, pvalue < spvalue)

