### Type I error ratio, average statistical power, and statistical powers of the four QTNs in three simulation scenarios.
#' @param gen represents the genotypes' data
#' @param matrix_final represents results multi-trait multi-locus GWAS method (MTML)
#' @param N is the execution/loop count
#' @param spvalue is the P-value threshold for significance
#' @param snp is the number of genotypes
#' @param r is the simulation scenarios
#' @param trait_k is the number of traits

statistic <- function(gen, matrix_final, N, spvalue, snp, r, trait_k) {
  genRaw <- as.matrix(gen[, 1:3])
  position <- genRaw[, 3]
  position <- as.numeric(position)
  QTN <- NULL
  QTN1 <- which(position >= (11298364 - 2000) & position <= (11298364 + 2000)) # 278
  QTN2 <- which(position >= (5066968 - 2000) & position <= (5066968 + 2000)) # 2054
  QTN3 <- which(position >= (5134228 - 2000) & position <= (5134228 + 2000)) # 2143
  QTN4 <- which(position >= (6137189 - 2000) & position <= (6137189 + 2000)) # 3716
  QTN <- rbind(QTN1, QTN2, QTN3, QTN4)

  positive_num <- vector()
  true_positive <- vector()
  true_negative <- vector()
  false_positive <- vector()
  false_negative <- vector()
  fp_ratio <- vector()
  fn_ratio <- vector()
  QTN_true_positive <- vector()
  QTN_false_negative <- vector()
  QTN_fn_ratio <- vector()

  a <- matrix(0, N, nrow(QTN))
  for (j in 1:nrow(QTN)) {
    for (i in 1:N) {
      last <- matrix_final[, 4:ncol(matrix_final)]
      last <- as.data.frame(last)
      positive_num[i] <- sum(last[, i] < spvalue)
      for (t in 1:nrow(last)) {
        if (last[t, i] < spvalue) {
          last[t, i] <- matrix_final[t, "Id"]
        }
      }
      a[i, j] <- length(intersect(last[, i], QTN[j, ]))
      if ((a[i, j] != 1) & (a[i, j] != 0)) a[i, j] <- 1
      true_positive[i] <- sum(a[i, ])
      QTN_true_positive[j] <- sum(a[, j])
      false_positive[i] <- positive_num[i] - true_positive[i]
      false_negative[i] <- nrow(QTN) - true_positive[i]
      QTN_false_negative[j] <- N - QTN_true_positive[j]
      true_negative[i] <- snp - true_positive[i] - false_positive[i] - false_negative[i]
      fp_ratio[i] <- false_positive[i] / (false_positive[i] + true_negative[i])
      fn_ratio[i] <- false_negative[i] / (false_negative[i] + true_positive[i])
      QTN_fn_ratio[j] <- QTN_false_negative[j] / (QTN_false_negative[j] + QTN_true_positive[j])
    }
  }
  fp_ratio <- mean(fp_ratio)
  fn_ratio <- mean(fn_ratio)
  power <- 1 - fn_ratio
  QTN_power <- rep(1, nrow(QTN)) - QTN_fn_ratio
  perform <- matrix(c(fp_ratio, power, QTN_power), 1, nrow(QTN) + 2)
  rownames(perform) <- paste("Y", r, "_", trait_k, sep = "")
  colnames(perform) <- c("Type_I_error", "Average_power", "QTN1_power", "QTN2_power", "QTN3_power", "QTN4_power")
  return(perform)
}
