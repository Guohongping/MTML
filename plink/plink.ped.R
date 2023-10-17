plink.ped <- function(n,phen,gene,p,q){
  ped1 <- matrix("NA",n,1)
  colnames(ped1) <- "FID"
  ped2 <- matrix("NA",n,1)
  colnames(ped2) <- "IID"
  ped3 <- matrix(0,n,1)
  colnames(ped3) <- "PID"
  ped4 <- matrix(0,n,1)
  colnames(ped4) <- "MID"
  ped5 <- matrix(0,n,1)
  colnames(ped5) <- "SEX"
  ped6 <- matrix(phen[,1],n,1)
  colnames(ped6) <- "phen"
  ped7 <- matrix(t(gene),n,q)
  colnames(ped7) <- matrix(1:q,q,1)
  ped <- cbind(ped1,ped2,ped3,ped4,ped5,ped6,ped7)
}