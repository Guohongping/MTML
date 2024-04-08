plink.ph <- function(n, phen, p) {
  ph3 <- phen
  colnames(ph3) <- matrix(1:p, p, 1)
  ph <- ph3
  return(ph)
}
