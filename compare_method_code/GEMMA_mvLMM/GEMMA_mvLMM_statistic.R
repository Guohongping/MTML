#***************************************** GEMMA: power and type I error *************************************

source("./code/statistic.R")
N <- 10
snp <- 4000
spvalue <- 0.05
gemma <- read.table("./output/gemma1_21.assoc.txt")
trait_info <- matrix(c(2, 5, 10, 10, 5, 2), 3, 2)
genotype <- read.csv("./data/simu_data/snps.csv", header = T)
gen <- as.matrix(genotype[, 1:3])

GEMMA_perform <- NULL
for (k in 1:nrow(trait_info)) {
  for (r in 1:3) {
    trait_k <- trait_info[k, 1]
    mvlmm <- matrix(gemma[, 2], nrow(gemma), 1)
    for (i in 1:N) {
      fname <- paste("./output/gemma", r, "_", trait_k, i, ".assoc.txt", sep = "")
      mv_table <- read.table(fname)
      mv_table <- as.matrix(mv_table)
      mvlmm <- cbind(mvlmm, mv_table[, ncol(mv_table)])
    }
    result <- mvlmm[1, 2:ncol(mvlmm)]
    result_final <- matrix(as.numeric(mvlmm[-1, ]), nrow(mvlmm[-1, ]), N + 1)
    colnames(result_final) <- c("Id", result)
    genRaw <- as.matrix(gen[, 1:3])
    colnames(genRaw) <- c("Id", "Chromosome", "Position")
    xpp <- genRaw
    ypp <- result_final
    matrix_final <- merge(xpp, ypp, by = "Id", all.x = TRUE)
    matrix_final[is.na(matrix_final)] <- 9
    statis <- statistic(gen, matrix_final, N, spvalue, snp, r, trait_k)
    GEMMA_perform <- rbind(GEMMA_perform, statis)
  }
}
write.csv(GEMMA_perform, "./results/GEMMA_result.csv")

#*************************************************** mvLMM: power and type I error **************************************
snp <- 4000
spvalue <- 0.05
mv <- read.table("./output/mvlmm1_21.assoc.txt")
trait_info <- matrix(c(2, 5, 10, 10, 5, 2), 3, 2)
genotype <- read.csv("./data/simu_data/snps.csv", header = T)
gen <- as.matrix(genotype[, 1:3])

mvLMM_perform <- NULL
for (k in 1:nrow(trait_info)) {
  for (r in 1:3) {
    trait_k <- trait_info[k, 1]
    N <- trait_info[k, 2]
    mvlmm <- matrix(mv[, 2], nrow(mv), 1)
    for (i in 1:N) {
      fname <- paste("./output/mvlmm", r, "_", trait_k, i, ".assoc.txt", sep = "")
      mv_table <- read.table(fname)
      mv_table <- as.matrix(mv_table)
      mvlmm <- cbind(mvlmm, mv_table[, ncol(mv_table)])
    }
    result <- mvlmm[1, 2:ncol(mvlmm)]
    result_final <- matrix(as.numeric(mvlmm[-1, ]), nrow(mvlmm[-1, ]), N + 1)
    colnames(result_final) <- c("Id", result)
    genRaw <- as.matrix(gen[, 1:3])
    colnames(genRaw) <- c("Id", "Chromosome", "Position")
    xpp <- genRaw
    ypp <- result_final
    matrix_final <- merge(xpp, ypp, by = "Id", all.x = TRUE)
    matrix_final[is.na(matrix_final)] <- 9
    statis <- statistic(gen, matrix_final, N, spvalue, snp, r, trait_k)
    mvLMM_perform <- rbind(mvLMM_perform, statis)
  }
}

write.csv(mvLMM_perform, "./results/mvLMM_result.csv")
