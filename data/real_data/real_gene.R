library(dplyr)
library(MASS)
library(readxl)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(bigmemory)
library(tibble)

kk <- NULL
psmatrix <- NULL
svpal <- 0.05
CLO <- 8
snp <- 215947
spvalue <- 0.05
trait_k <- 3
N <- 1

genotype <- read.csv("./data/real_data/gen_deal.csv", header = T)
gen <- cbind(as.matrix(genotype[, 1:3]), as.matrix(genotype[, 4:202]))
genRaw <- as.matrix(gen[, 1:3])
phenotype <- read.csv("./data/real_data/trait_f3.csv", header = T)
phenotype <- as.matrix(phenotype[, 4:ncol(phenotype)])
gen_info <- read.csv("./data/real_data/gen_info.csv")
gen_info <- as.matrix(gen_info)
### add population structure
ph_vec_mat <- function(ph) {
  ph <- as.matrix(ph)
  ph_matrix <- matrix(ph, 199, trait_k)
  ph_matrix <- as.data.frame(ph_matrix)
  colnames(ph_matrix) <- paste("ph_", 1:trait_k)
  return(ph_matrix = ph_matrix)
}
phen <- ph_vec_mat(phenotype)
phen <- as.matrix(phen)
str1 <- read.table("./data/real_data/real_ph_f3.3.Q", header = F)
ps <- as.matrix(str1)

if (is.null(ps) == FALSE) {
  ps1 <- cbind(matrix(1, nrow = nrow(phen)), ps)
  vhat <- solve(crossprod(ps1, ps1)) %*% crossprod(ps1, phen)
  vhat1 <- vhat[-1, ]
  phen_ph <- phen - ps %*% vhat1
} else {
  phen_ph <- phen
}
phen_ph <- as.matrix(phen_ph)

#******************************************* MTML*******************************************************

source("./code/single_trait.R")
source("./code/cauchyComb.R")
source("./code/ridge.R")
source("./code/MTML.R")

matrix_final <- MTML(gen, phen_ph, kk, psmatrix, N, svpal, CLO, trait_k, snp)
matrix_final <- na.omit(matrix_final)
colnames(matrix_final) <- c("Id", "Chromosome", "Position", "pvalue")
matrix_final <- data.frame(matrix_final)
matrix_final[, "pvalue"] <- as.numeric(matrix_final[, "pvalue"])
trait_final <- filter(matrix_final, pvalue < spvalue)
trait_final <- as.data.frame(trait_final)
mtml <- rownames_to_column(trait_final, var = "order")

mtml_detect_gene <- NULL
for (k in 1:nrow(mtml)) {
  for (m in 1:nrow(gen_info)) {
    if ((mtml[k, 3] == gen_info[m, 1]) & (as.numeric(mtml[k, 4]) >= as.numeric(gen_info[m, 2]) - 20000) & (as.numeric(mtml[k, 4]) <= as.numeric(gen_info[m, 3]) + 20000)) {
      mtml_gene <- c(mtml[k, 1], mtml[k, 2], mtml[k, 3], mtml[k, 4], gen_info[m, 4], mtml[k, 5])
      mtml_gene <- as.matrix(mtml_gene, 1, 6)
      mtml_detect_gene <- cbind(mtml_detect_gene, mtml_gene)
    }
  }
}
mtml_detect_gene <- t(mtml_detect_gene)
colnames(mtml_detect_gene) <- c("order", "ID", "Chromosome", "Position", "gene", "pvalue")
mtml_detect_gene_final <- mtml_detect_gene[!duplicated(mtml_detect_gene[, 2]), ]
mtml_detect_gene_final <- mtml_detect_gene_final[!duplicated(mtml_detect_gene_final[, 5]), ]
mtml_gene_number <- nrow(mtml_detect_gene_final)
write.csv(mtml_detect_gene_final, "./results/MTML_real_result.csv")


#******************************************* FASTmrEMMA*******************************************************

source("./compare_method_code/FASTmrEMMA/myFASTmrEMMA.R")

phe <- as.matrix(phen_ph[, 1])
result <- myFASTmrEMMA(gen = gen, phe = phe, genRaw = genRaw, kk = NULL, psmatrix = NULL, svpal = 0.05, svmlod = 3, CLO = 16)
result1 <- as.data.frame(result)
trait1 <- rownames_to_column(result1, var = "order")
tt <- c(rep("trait1", times = nrow(trait1)))
trait1 <- cbind(tt, trait1)
phe <- as.matrix(phen_ph[, 2])
result <- myFASTmrEMMA(gen = gen, phe = phe, genRaw = genRaw, kk = NULL, psmatrix = NULL, svpal = 0.05, svmlod = 3, CLO = 16)
result2 <- as.matrix(as.data.frame(result))
result2 <- as.data.frame(t(result2))
if (nrow(result2) == 1) {
  trait2 <- NULL
}else {
  result2 <- as.data.frame(result)
  trait2 <- rownames_to_column(result2, var = "order")
  tt <- c(rep("trait2", times = nrow(trait2)))
  trait2 <- cbind(tt, trait2)
}
phe <- as.matrix(phen_ph[, 3])
result <- myFASTmrEMMA(gen = gen, phe = phe, genRaw = genRaw, kk = NULL, psmatrix = NULL, svpal = 0.05, svmlod = 3, CLO = 16)
result3 <- as.data.frame(result)
trait3 <- rownames_to_column(result3, var = "order")
tt <- c(rep("trait3", times = nrow(trait3)))
trait3 <- cbind(tt, trait3)
fast <- rbind(trait1, trait2, trait3)
fast <- as.matrix(fast)
fast_detect_gene <- NULL
for (k in 1:nrow(fast)) {
  for (m in 1:nrow(gen_info)) {
    if ((fast[k, 5] == gen_info[m, 1]) & (as.numeric(fast[k, 3]) > as.numeric(gen_info[m, 2]) - 20000) & (as.numeric(fast[k, 3]) < as.numeric(gen_info[m, 3]) + 20000)) {
      fast_gene <- c(as.character(fast[k, 1]), fast[k, 2], fast[k, 4], fast[k, 5], fast[k, 3], gen_info[m, 4])
      fast_gene <- as.matrix(fast_gene, nrow = 1, ncol = 6)
      fast_detect_gene <- cbind(fast_detect_gene, fast_gene)
    }
  }
}
fast_detect_gene <- t(fast_detect_gene)
colnames(fast_detect_gene) <- c("trait_name", "order", "ID", "Chromosome", "Position", "gene")
fast_detect_gene_final <- fast_detect_gene[!duplicated(fast_detect_gene[, 3]), ]
fast_detect_gene_final <- fast_detect_gene_final[!duplicated(fast_detect_gene_final[, 6]), ]

fast_gene <- NULL
for (i in 1:trait_k) {
  fname <- paste("fast_trait", i, sep = "")
  name <- paste("trait", i, sep = "")
  fname <- as.character(fname)
  name <- as.character(name)
  if (name %in% fast_detect_gene_final[, 1] == TRUE) {
    fname <- sum(fast_detect_gene_final[, 1] == name)
  } else {
    fname <- 0
  }
  fast_gene <- cbind(fast_gene, fname)
}

fast_gene_number <- fast_gene
colnames(fast_gene_number) <- c("trait1", "trait2", "trait3")
write.csv(fast_detect_gene_final, "./results/FASTmrEMMA_real_result.csv")

#******************************************* GEMMA**************************************************************

spvalue <- 0.00000023
gemma1 <- read.table("./output/gemmaf31.assoc.txt", header = T)
gemma <- cbind(gemma1[, 2], gemma1[, ncol(gemma1)])
colnames(gemma) <- c("Id", "pvalue")
gemma <- data.frame(gemma)
gemma[, "pvalue"] <- as.numeric(gemma[, "pvalue"])
trait1 <- filter(gemma, pvalue < spvalue)
trait1 <- rownames_to_column(trait1, var = "order")
tt <- c(rep("trait1", times = nrow(trait1)))
trait1 <- cbind(tt, trait1)
trait1_final <- merge(trait1, genRaw, by = "Id", all.x = TRUE)

gemma2 <- read.table("./output/gemmaf32.assoc.txt", header = T)
gemma <- cbind(gemma2[, 2], gemma2[, ncol(gemma1)])
colnames(gemma) <- c("Id", "pvalue")
gemma <- data.frame(gemma)
gemma[, "pvalue"] <- as.numeric(gemma[, "pvalue"])
trait2 <- filter(gemma, pvalue < spvalue)
trait2 <- rownames_to_column(trait2, var = "order")
tt <- c(rep("trait2", times = nrow(trait2)))
trait2 <- cbind(tt, trait2)
trait2_final <- merge(trait2, genRaw, by = "Id", all.x = TRUE)

gemma3 <- read.table("./output/gemmaf33.assoc.txt", header = T)
gemma <- cbind(gemma3[, 2], gemma3[, ncol(gemma1)])
colnames(gemma) <- c("Id", "pvalue")
gemma <- data.frame(gemma)
gemma[, "pvalue"] <- as.numeric(gemma[, "pvalue"])
trait3 <- filter(gemma, pvalue < spvalue)
trait3 <- rownames_to_column(trait3, var = "order")
tt <- c(rep("trait3", times = nrow(trait3)))
trait3 <- cbind(tt, trait3)
trait3_final <- merge(trait3, genRaw, by = "Id", all.x = TRUE)

gemma <- rbind(trait1_final, trait2_final, trait3_final)
gemma <- cbind(gemma[, 2:3], gemma[, 1], gemma[, 5:6], gemma[, 4])
gemma_detect_gene <- NULL
for (k in 1:nrow(gemma)) {
  for (m in 1:nrow(gen_info)) {
    if ((gemma[k, 4] == gen_info[m, 1]) & (as.numeric(gemma[k, 5]) > as.numeric(gen_info[m, 2]) - 20000) & (as.numeric(gemma[k, 5]) < as.numeric(gen_info[m, 3]) + 20000)) {
      gemma_gene <- c(as.character(gemma[k, 1]), gemma[k, 2], gemma[k, 3], gemma[k, 4], gemma[k, 5], gen_info[m, 4])
      gemma_gene <- as.matrix(gemma_gene, nrow = 1, ncol = 6)
      gemma_detect_gene <- cbind(gemma_detect_gene, gemma_gene)
    }
  }
}
gemma_detect_gene <- t(gemma_detect_gene)
colnames(gemma_detect_gene) <- c("trait_name", "order", "ID", "Chromosome", "Position", "gene")
gemma_detect_gene_final <- gemma_detect_gene[!duplicated(gemma_detect_gene[, 3]), ]
gemma_detect_gene_final <- gemma_detect_gene_final[!duplicated(gemma_detect_gene_final[, 6]), ]

gemma_gene <- NULL
for (i in 1:trait_k) {
  fname <- paste("gemma_trait", i, sep = "")
  name <- paste("trait", i, sep = "")
  fname <- as.character(fname)
  name <- as.character(name)
  if (name %in% gemma_detect_gene_final[, 1] == TRUE) {
    fname <- sum(gemma_detect_gene_final[, 1] == name)
  } else {
    fname <- 0
  }
  gemma_gene <- cbind(gemma_gene, fname)
}

gemma_gene_number <- gemma_gene
colnames(gemma_gene_number) <- c("trait1", "trait2", "trait3")
write.csv(gemma_detect_gene_final, "./results/GEMMA_real_result.csv")


#******************************************* mvLMM**************************************************************
spvalue <- 0.00000023
mvlmm <- read.table("./output/mvlmmf3.assoc.txt", header = T)
mvlmm <- cbind(mvlmm[, 2], mvlmm[, ncol(mvlmm)])
colnames(mvlmm) <- c("Id", "pvalue")
mvlmm <- data.frame(mvlmm)
mvlmm[, "pvalue"] <- as.numeric(mvlmm[, "pvalue"])
mvlmm <- filter(mvlmm, pvalue < spvalue)
mvlmm_gene_number <- sum(as.numeric(mvlmm[, 2] != 0))
if (mvlmm_gene_number == 0) {
  mvlmm_detect_gene <- c(0, 0, 0, 0, 0)
  mvlmm_detect_gene_final <- t(mvlmm_detect_gene)
  colnames(mvlmm_detect_gene_final) <- c("order", "ID", "Chromosome", "Position", "gene")
} else {
  mvlmm_detect_gene <- NULL
  for (k in 1:nrow(mvlmm)) {
    for (m in 1:nrow(gen_info)) {
      if ((mvlmm[k, 3] == gen_info[m, 1]) & (as.numeric(mvlmm[k, 4]) >= as.numeric(gen_info[m, 2]) - 20000) & (as.numeric(mvlmm[k, 4]) <= as.numeric(gen_info[m, 3]) + 20000)) {
        mvlmm_gene <- c(mvlmm[k, 1], mtml[k, 2], mvlmm[k, 3], mvlmm[k, 4], gen_info[m, 4])
        mvlmm_gene <- as.matrix(mvlmm_gene, 1, 5)
        mvlmm_detect_gene <- cbind(mvlmm_detect_gene, mvlmm_gene)
      }
    }
  }
  mvlmm_detect_gene <- t(mvlmm_detect_gene)
  colnames(mvlmm_detect_gene) <- c("order", "ID", "Chromosome", "Position", "gene")
  mvlmm_detect_gene_final <- mvlmm_detect_gene[!duplicated(mvlmm_detect_gene[, 2]), ]
  mvlmm_detect_gene_final <- mvlmm_detect_gene_final[!duplicated(mvlmm_detect_gene_final[, 5]), ]
}
write.csv(mvlmm_detect_gene_final, "./results/mvLMM_real_result.csv")

#******************************************** TABLE1 **********************************************************
detect_real_gene_number <- rbind(mtml_gene_number, gemma_gene_number, fast_gene_number, mvlmm_gene_number)
colnames(detect_real_gene_number) <- c("trait1", "trait2", "trait3")
rownames(detect_real_gene_number) <- c("MTML", "GEMMA", "FASTmrEMMA", "mvLMM")
write.csv(detect_real_gene_number, "./results/Table1.csv")
