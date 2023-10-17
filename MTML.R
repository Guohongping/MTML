 ### MTML: An efficient multi-trait multi-locus GWAS method based on Cauchy combination test
 #' three Steps of multi-trait multi-locus GWAS method (MTML)
 #' single-trait: matrix transformation for single-trait decomposition analyses
 #' multi-locus: dimension reduction and deshrinking ridge regression for multi-locus analyses
 #' multi-trait: cauchy combination test for aggregating the individual p-values of multiple traits
 
 #' @param gen represents the genotypes' data
 #' @param phenotype represents the traits' data
 #' @param kk is the initial kinship matrix
 #' @param psmatrix is the initial P-value
 #' @param N is the execution/loop count
 #' @param svpal is the P-value threshold for significance
 #' @param Genformat is the gene mode (including 0 and 1)
 #' @param Likelihood is the likelihood function type, "REML" represents restrictive maximum likelihood; "ML" represents maximum likelihood
 #' @param CLO is the the number of computer cores
 #' @param snp is the number of genotypes
 #' @param trait_k is the number of traits
 
 MTML <- function(gen,phenotype,kk,psmatrix,N,svpal,Genformat,Likelihood,CLO,trait_k,snp){
 
    matrix_final <- gen[,1:3]
    for (t in 1:N){
	
    Y <- as.matrix(phenotype[,((trait_k)*(t-1)+1):((trait_k)*t)])
	  final <- NULL
    xp <- gen[,1:3]
    for (i in 1:trait_k){
      phe <- as.matrix(Y[,i]) 
      result <- single_trait(gen=gen,phe=phe,kk=NULL,psmatrix=NULL,svpal,Genformat,Likelihood="REML",CLO)
      result <- as.data.frame(result)
      result_i <- result[,4]
      result_i <- 10^(-abs(as.numeric(result_i)))
      matrix_p <- cbind(gen[,1:3],result_i)
      colnames(matrix_p) <- c("ID","Chromosome","Position","sspval")
	  matrix_p <- as.data.frame(matrix_p)
	  ssp <- filter(matrix_p,sspval < 0.1)
	     
      xpp <- ssp[1:3]
      ypp <- gen[,3:ncol(gen)]
      gen0 <- merge(xpp,ypp,by = "Position",all.x = TRUE)
      genraw <- gen0[,1:3]
      gen0<-as.matrix(gen0[,4:ncol(gen0)])
      gen0 <- t(gen0)
      m<-ncol(gen0)
      n<-nrow(gen0)
      kk0<-matrix(0,n,n)
      for(k in 1:m){
        z <- as.matrix(gen0[,k])
        kk0 <- kk0+z%*%t(z)
      }
      row.names(kk0)<-NULL
      rm(z)
      x1<-rep(1,times=n)
      x<-as.matrix(x1)
      z<-gen0
      qq<-eigen(kk0,symmetric=T)
      
      y <- phe
      fit<-ridge(x=x,y=y,z=z,qq=qq)
      vg<-fit[[1]]$vg
      ve<-fit[[1]]$ve
      lambda<-fit[[1]]$lambda
      gamma<-fit[[3]]$gamma
      vgk1<-fit[[3]]$vgk1
      vgk2<-fit[[3]]$vgk2
      df<-fit[[3]]$df
      wald<-gamma^2/vgk1
      w.wald<-gamma^2/vgk2
      d.gamma<-gamma/df
      d.wald<-wald/df
      test<-data.frame(vg,ve,lambda,gamma,wald,df,d.gamma,d.wald,w.wald)
      names(test)<-c("vg","ve","lambda","gamma","ORR.wald","df","DRR.gamma","DRR.wald","HRR.wald")
      ORR.wald<-test$ORR.wald
      DRR.wald<-test$DRR.wald
      HRR.wald<-test$HRR.wald
      
      ORR.p<-pchisq(ORR.wald,1,lower.tail=F)
      DRR.p<-pchisq(DRR.wald,1,lower.tail=F)
      HRR.p<-pchisq(HRR.wald,1,lower.tail=F)
	  
	  matrix_s <- cbind(genraw,DRR.p)
	  slpval <- as.matrix(cbind(matrix_s[,1],matrix_s[,ncol(matrix_s)]))
	  colnames(slpval) <- c("Position","slpval")
	    
      yp <- slpval
      slp <- merge(xp,yp,by = "Position")
	  xp <- slp
    }
    
	final <- as.matrix(xp[,4:ncol(xp)])
  
    pval <- cauchyComb(final,weig = NULL,nr = nrow(final) ,nc = ncol(final))
    final_p1 <- cbind(xp[,1:3],pval)
    final_p1 <- data.frame(final_p1)
    final_p1[,"pval"] <- as.numeric(final_p1[,"pval"])
    mt_pval <- filter(final_p1,pval < svpal)
    mt_pval <- as.matrix(cbind(mt_pval[,1],mt_pval[,4]))
    colnames(mt_pval) <- c("Position","pvalue")


    xpp <- gen[,1:3]
    ypp <- mt_pval
    final_p2 <- merge(xpp,ypp,by = "Position",all.x = TRUE)
    final_p3 <- final_p2[order(final_p2$Id),]
    matrix_final <- cbind(matrix_final,final_p3[,4])
    }

	return (matrix_final)
}