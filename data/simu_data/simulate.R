### Three Monte Carlo simulations
#' @param trait_k is the number of traits
#' @param N is the execution/loop count
#' @param X represents the genotypes' data for four supposed pleiotropic QTNs 
#' @param genetic_effect1 is the genetic effect value for four supposed pleiotropic QTNs in simulation 1 and simulation 2
#' @param genetic_effect2 is the genetic effect value for four supposed pleiotropic QTNs in simulation 3
#' @param trait_info is the information for traits


simulate <- function(trait_k,N,X,genetic_effect1,genetic_effect2){
  Y1 <- NULL
  Y2 <- NULL
  Y3 <- NULL
  mu <- matrix(rep(1,nrow(X))) %*% t(rep(10,trait_k))  
  
  b1 <- matrix(rep(genetic_effect1,trait_k),nrow=4,ncol=trait_k)
  b2 <- matrix(rep(genetic_effect2,trait_k),nrow=4,ncol=trait_k)   
  
  V_e_0 <- diag(10,trait_k)
  V_e <- diag(10,trait_k)    
  V_e[V_e == 0] <- 5 
  V_g <- diag(2,trait_k)     
  V_g[V_g == 0] <- 1

  for (i in 1:N){
    set.seed(i)
	E_0 <- rmvnorm(nrow(X),mean=rep(0,trait_k),sigma=V_e_0)
    E <- rmvnorm(nrow(X),mean=rep(0,trait_k),sigma=V_e)
	G <- rmvnorm(nrow(X),mean=rep(0,trait_k),sigma=V_g)
    y1 <- mu + X%*%b1 + E_0
	y2 <- mu + X%*%b1 + E
	y3 <- mu + X%*%b2 + G + E
	
	if (trait_k == 2){
	name1_1 <- paste("trait1-",i,sep="")
	name1_2 <- paste("trait2-",i,sep="")
	colnames(y1) <- c(name1_1,name1_2)
	colnames(y2) <- c(name1_1,name1_2)
	colnames(y3) <- c(name1_1,name1_2)
	}else if (trait_k ==5){
	name3_1 <- paste("trait1-",i,sep="")
	name3_2 <- paste("trait2-",i,sep="")
	name3_3 <- paste("trait3-",i,sep="")
	name3_4 <- paste("trait4-",i,sep="")
	name3_5 <- paste("trait5-",i,sep="")
	colnames(y1) <- c(name3_1,name3_2,name3_3,name3_4,name3_5)
	colnames(y2) <- c(name3_1,name3_2,name3_3,name3_4,name3_5)
	colnames(y3) <- c(name3_1,name3_2,name3_3,name3_4,name3_5)
	}else if (trait_k ==10){
	name4_1 <- paste("trait1-",i,sep="")
	name4_2 <- paste("trait2-",i,sep="")
	name4_3 <- paste("trait3-",i,sep="")
	name4_4 <- paste("trait4-",i,sep="")
	name4_5 <- paste("trait5-",i,sep="")
	name4_6 <- paste("trait6-",i,sep="")
	name4_7 <- paste("trait7-",i,sep="")
	name4_8 <- paste("trait8-",i,sep="")
	name4_9 <- paste("trait9-",i,sep="")
	name4_10 <- paste("trait10-",i,sep="")
	colnames(y1) <- c(name4_1,name4_2,name4_3,name4_4,name4_5,name4_6,name4_7,name4_8,name4_9,name4_10)
	colnames(y2) <- c(name4_1,name4_2,name4_3,name4_4,name4_5,name4_6,name4_7,name4_8,name4_9,name4_10)
	colnames(y3) <- c(name4_1,name4_2,name4_3,name4_4,name4_5,name4_6,name4_7,name4_8,name4_9,name4_10)
	}
    Y1 <- cbind(Y1,y1)
	Y2 <- cbind(Y2,y2)
	Y3 <- cbind(Y3,y3)
  }
  Y <- list(Y1,Y2,Y3)
  return(result = Y)
}




simu_data <- function(QTN,trait_info,genetic_effect1,genetic_effect2){
  X <- as.matrix(QTN[,4:ncol(QTN)])
  X <- t(X)  
  for (j in 1:nrow(trait_info)){
    trait_k <- trait_info[j,1]
    N <- trait_info[j,2]
    simu_yy <- simulate(trait_k,N,X,genetic_effect1,genetic_effect2)
    file1 <- simu_yy[1]
    file2 <- simu_yy[2]
    file3 <- simu_yy[3]
    rname <- seq(1,nrow(X),1)
    file1 <- as.data.frame(file1)
    file2 <- as.data.frame(file2)
    file3 <- as.data.frame(file3)
    rownames(file1) <- paste("Ind",rname,sep="-")
    rownames(file2) <- paste("Ind",rname,sep="-")
    rownames(file3) <- paste("Ind",rname,sep="-")
    fname1 <- paste("./MTML_code_data/data/simu_data/simu_phenotype/Y1_",trait_k,".csv",sep="")
    fname2 <- paste("./MTML_code_data/data/simu_data/simu_phenotype/Y2_",trait_k,".csv",sep="")
    fname3 <- paste("./MTML_code_data/data/simu_data/simu_phenotype/Y3_",trait_k,".csv",sep="")
    write.csv(file1,fname1)
    write.csv(file2,fname2)
    write.csv(file3,fname3)
  }
}  