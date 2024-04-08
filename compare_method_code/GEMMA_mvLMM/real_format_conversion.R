source("./compare_method_code/GEMMA_mvLMM/plink.map.R")
source("./compare_method_code/GEMMA_mvLMM/plink.ped.R")
source("./compare_method_code/GEMMA_mvLMM/plink.ph.R")

genotype <- read.csv("./data/real_data/gen_deal.csv", header = T)
genotype <- as.matrix(genotype[, 4:202])
genotype[genotype == 0] <- "11"
genotype[genotype == 1] <- "12"
genotype[genotype == 2] <- "22"
gene <- genotype
gene <- t(gene)
phenotype <- read.csv("./data/real_data/trait_f3.csv", header = T)
phenotype <- as.matrix(phenotype[, 4:ncol(phenotype)])
trait_k <- 3
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
phen <- phen_ph

n <- nrow(gene)
q <- ncol(gene)
p <- ncol(phen)
mapph <- plink.map(q)
write.table(mapph, "./data/real_data/trait_f3.map", row.names = F, col.names = F, quote = F)
peddh <- plink.ped(n, phen, gene, p, q)
write.table(peddh, "./data/real_data/trait_f3.ped", row.names = F, col.names = F, quote = F)
phh <- plink.ph(n, phen, p)
write.table(phh, "./data/real_data/trait_f3.txt", row.names = F, col.names = F, quote = F)
