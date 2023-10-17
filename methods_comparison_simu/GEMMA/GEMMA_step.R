
R
source("plink.map.R")
source("plink.ped.R")
source("plink.ph.R")

genotype <- read.csv("snps.csv",header=T)
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
    fname <- paste("Y",r,"_",trait_k,".csv",sep="")
    phenotype <- read.csv(fname,header = T)
    phenotype <- as.matrix(phenotype[,-1])
	phen <- phenotype
	
	n <- nrow(gene) 
    q <- ncol(gene) 
    p <- ncol(phen)
    mapph <- plink.map(q)
    peddh <- plink.ped(n,phen,gene,p,q)
	phh <- plink.ph(n,phen,p)
	
	map_fname <- paste("y",r,"_",trait_k,".map",sep="")
	ped_fname <- paste("y",r,"_",trait_k,".ped",sep="")
	ph_fname <-  paste("y",r,"_",trait_k,".txt",sep="")
	write.table(mapph,map_fname,row.names=F,col.names=F,quote=F)
	write.table(peddh,ped_fname,row.names=F,col.names=F,quote=F)
	write.table(pph,ph_fname,row.names=F,col.names=F,quote=F)
  }
}




#****************************************** Plink bed et al ************************************************
./plink --file y1_2 --make-bed --out y1_2

#******************************************** Gemma related matrix************************************************
./gemma-0.98.5-linux-static-AMD64 -bfile y1_2 -p ph1_2.txt -gk 1 -o related1_2

# #******************************************* Gemma ************************************************
for ((i=1;i<=10;i++))
do
  ./gemma-0.98.5-linux-static-AMD64 -bfile y1_2 -p ph1_2.txt  -k related1_2.cXX.txt -lmm 1 -n $i -o gemma1_2$i
  
done

#***************************************** power and type I error ******************************************************************
R
source("statistic.R")
N <- 10
snp <- 4000
spvalue <- 0.05
gemma <- read.table("gemma1_21.assoc.txt")
trait_info <- matrix(c(2,5,10,10,4,2),3,2)
genotype <- read.csv("snps.csv",header=T)
gen <- as.matrix(genotype[,1:3])
	
GEMMA_perform <- NULL
for (k in 1:nrow(trait_info)){
  for (r in 1:3){
    trait_k <- trait_info[k,1]
	mvlmm <- matrix(gemma[,2],nrow(gemma),1)
    for (i in 1:N){
      fname <- paste("gemma",r,"_",trait_k,i,".assoc.txt",sep="")
      mv_table <- read.table(fname)
      mv_table <- as.matrix(mv_table)
      mvlmm <- cbind(mvlmm,mv_table[,ncol(mv_table)])
    }
    result <- mvlmm[1,2:ncol(mvlmm)]
    result_final <- matrix(as.numeric(mvlmm[-1,]),nrow(mvlmm[-1,]),N+1)
    colnames(result_final) <- c("Id",result)
    genRaw <- as.matrix(gen[,1:3])
    colnames(genRaw) <- c("Id","Chromosome","Position")
    xpp <- genRaw
    ypp <- result_final
    matrix_final <- merge(xpp,ypp,by="Id",all.x = TRUE)
    matrix_final[is.na(matrix_final)] = 9
    statis <- statistic(gen,matrix_final,N,spvalue,snp,r,trait_k)
    GEMMA_perform <- rbind(GEMMA_perform,statis)
  }
}




