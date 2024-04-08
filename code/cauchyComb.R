### Cauchy combination test for aggregating the individual p-values of multiple traits
#' @param weig represents the weight coefficient
#' @param final is the P-value obtained from single trait multi loci analysis
#' @param nr is the number of genotypes
#' @param nc is the number of traits

cauchyComb <- function(final, weig = NULL, nr, nc) {
  if (is.null(weig) == TRUE) {
    wei <- ncol(final)
    weig <- 1 / wei
  }
  sta0 <- rowSums(matrix((tan((0.5 - as.numeric(final)) * pi) * weig), nr, nc))
  pval <- 0.5 - atan(sta0) / pi
  return(pval)
}
