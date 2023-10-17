
# #******************************************* mvLMM ************************************************
#d=2
for ((i=1;i<=5;i++))
do
  ./gemma-0.98.5-linux-static-AMD64 -bfile y1_2 -p ph1_2.txt  -k related1_2.cXX.txt  -lmm 1 -n $[2$i-1] $[2$i] -o mvlmm1_2$i  
done



#d=5
for ((i=1;i<=200;i++))
do
  ./gemma-0.98.5-linux-static-AMD64 -bfile y1_5 -p ph1_5.txt  -k related1_5.cXX.txt  -lmm 1 -n $[5$i-4] $[5$i-3] $[5$i-2] $[5$i-1] $[5$i] -o mvlmm1_5$i  
done


#d=10
for ((i=1;i<=100;i++))
do
  ./gemma-0.98.5-linux-static-AMD64 -bfile y1_10 -p ph1_10.txt  -k related1_10.cXX.txt  -lmm 1 -n $[10$i-9] $[10$i-8] $[10$i-7] $[10$i-6] $[10$i-5] $[10$i-4] $[10$i-3] $[10$i-2] $[10$i-1] $[10$i] -o mvlmm1_10$i  
done


#***************************************** power and type I error ******************************************************************


R

snp <- 4000
spvalue <- 0.05
mv <- read.table("mvlmm1_21.assoc.txt")
trait_info <- matrix(c(2,5,10,10,4,2),3,2)
genotype <- read.csv("snps.csv",header=T)
gen <- as.matrix(genotype[,1:3])
	
mvLMM_perform <- NULL
for (k in 1:nrow(trait_info)){
  for (r in 1:3){
    trait_k <- trait_info[k,1]
	N <- trait_info[k,2]
	mvlmm <- matrix(mv[,2],nrow(mv),1)
    for (i in 1:N){
      fname <- paste("mvlmm",r,"_",trait_k,i,".assoc.txt",sep="")
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
    mvLMM_perform <- rbind(mvLMM_perform,statis)
  }
}


