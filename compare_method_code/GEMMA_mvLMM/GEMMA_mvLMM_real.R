
library(dplyr)
spvalue <- 0.00000023
gemma1 <- read.table("./output/gemmaf31.assoc.txt",header=T)
gemma<- cbind(gemma1[,2],gemma1[,ncol(gemma1)])
colnames(gemma)<- c("Id","pvalue")
gemma <- data.frame(gemma)
gemma[,"pvalue"] <- as.numeric(gemma[,"pvalue"])
trait1 <- filter(gemma, pvalue < spvalue) 


gemma2 <- read.table("./output/gemmaf32.assoc.txt",header=T)
gemma<- cbind(gemma2[,2],gemma2[,ncol(gemma1)])
colnames(gemma)<- c("Id","pvalue")
gemma <- data.frame(gemma)
gemma[,"pvalue"] <- as.numeric(gemma[,"pvalue"])
trait2 <- filter(gemma, pvalue < spvalue) 

gemma3 <- read.table("./output/gemmaf33.assoc.txt",header=T)
gemma<- cbind(gemma3[,2],gemma3[,ncol(gemma1)])
colnames(gemma)<- c("Id","pvalue")
gemma <- data.frame(gemma)
gemma[,"pvalue"] <- as.numeric(gemma[,"pvalue"])
trait3 <- filter(gemma, pvalue < spvalue) 
trait_result <- rbind(trait1,trait2,trait3)
write.csv(trait_result,"./results/GEMMA_mvLMM_real_result.csv")

