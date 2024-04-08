plink.map <- function(q) {
  map1 <- matrix(0, q, 1)
  colnames(map1) <- "chr"
  map2 <- matrix(1:nrow(map1), nrow(map1), 1)
  colnames(map2) <- "marker name"
  map3 <- matrix("NA", q, 1)
  colnames(map3) <- "genetic distance"
  map4 <- matrix(1:nrow(map1), nrow(map1), 1)
  colnames(map4) <- "position"
  map <- cbind(map1, map2, map3, map4)
}
