source("./compare_method_code/GEMMA_mvLMM/plink.map.R")
source("./compare_method_code/GEMMA_mvLMM/plink.ped.R")
source("./compare_method_code/GEMMA_mvLMM/plink.ph.R")

genotype <- read.csv("./data/simu_data/snps.csv",header=T)
genotype <- as.matrix(genotype[,4:202])
genotype[genotype==0] <- "11"
 genotype[genotype==1] <- "12"
 genotype[genotype==2] <- "22"
 gene <- genotype
 gene <- t(gene)
trait_info <- matrix(c(2,5,10,10,4,2),3,2)

for (k in 1:nrow(trait_info)){
  for (r in 1:3){
    trait_k <- trait_info[k,1]
    N <- trait_info[k,2]
    fname <- paste("./data/simu_data/simu_phenotype/Y",r,"_",trait_k,".csv",sep="")
    phenotype <- read.csv(fname,header = T)
    phenotype <- as.matrix(phenotype[,-1])
	phen <- phenotype
	
	n <- nrow(gene) 
    q <- ncol(gene) 
    p <- ncol(phen)
    mapph <- plink.map(q)
    peddh <- plink.ped(n,phen,gene,p,q)
	phh <- plink.ph(n,phen,p)
	
	map_fname <- paste("./compare_method_code/GEMMA_mvLMM/y",r,"_",trait_k,".map",sep="")
	ped_fname <- paste("./compare_method_code/GEMMA_mvLMM/y",r,"_",trait_k,".ped",sep="")
	ph_fname <-  paste("./compare_method_code/GEMMA_mvLMM/y",r,"_",trait_k,".txt",sep="")
	write.table(mapph,map_fname,row.names=F,col.names=F,quote=F)
	write.table(peddh,ped_fname,row.names=F,col.names=F,quote=F)
	write.table(phh,ph_fname,row.names=F,col.names=F,quote=F)
  }
}

