library(dplyr)
library(MASS)
library(readxl)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(bigmemory)

source("./compare_method_code/FASTmrEMMA/myFASTmrEMMA.R")
source("./code/statistic.R")

kk <- NULL
psmatrix <- NULL
Likelihood <- "REML"
svpal <- 0.05
svmlod <- 3
Genformat <- 1
CLO <- 8
N <- 10
snp <- 4000
spvalue <- 0.0002

genotype <- read.csv("./data/simu_data/snps.csv", header = T)
gen <- cbind(as.matrix(genotype[, 1:3]), as.matrix(genotype[, 4:202]))
genRaw <- as.matrix(gen[, 1:3])
trait_info <- matrix(c(2, 5, 10, 10, 4, 2), 3, 2)

FASTmrEMMA_perform <- NULL
for (k in 1:nrow(trait_info)) {
  for (r in 1:3) {
    trait_k <- trait_info[k, 1]
    fname <- paste("./data/simu_data/simu_phenotype/Y", r, "_", trait_k, ".csv", sep = "")
    phenotype <- read.csv(fname, header = T)
    phenotype <- as.matrix(phenotype[, -1])

    result5 <- NULL
    for (i in 1:N) {
      phe <- as.matrix(phenotype[, i])
      result <- myFASTmrEMMA(gen, phe, genRaw, kk, psmatrix, svpal, svmlod, CLO)
      result1 <- as.matrix(as.data.frame(result))
      result2 <- result1[, 6]
      result3 <- 10^(-abs(as.numeric(result2)))
      final <- cbind(result1[, 1], result3)
      colnames(final) <- c("Position", "pvalue")
      xpp <- gen[, 1:3]
      ypp <- final
      final_p2 <- merge(xpp, ypp, by = "Position", all.x = TRUE)
      final_p2 <- as.data.frame(final_p2)
      final_final <- final_p2[order(final_p2[, 2], decreasing = FALSE), ]
      result4 <- as.matrix(as.numeric(as.character(final_final[, 4])))
      result5 <- cbind(result5, result4)
    }
    fast_result <- cbind(genRaw[, 1:3], result5)
    matrix_final <- fast_result

    matrix_final[is.na(matrix_final)] <- 9
    statis <- statistic(gen, matrix_final, N, spvalue, snp, r, trait_k)
    FASTmrEMMA_perform <- rbind(FASTmrEMMA_perform, statis)
  }
}
write.csv(FASTmrEMMA_perform, "./results/FASTmrEMMA_result.csv")
