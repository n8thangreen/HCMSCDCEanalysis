A <- matrix(c(
  0.152233, 1.4594, 7.5275,
  0.146514, 1.3434, 7.6672,
  0.301198, 2.7718, 6.3117
), nrow = 3, byrow = TRUE)

B <- c(
  1899 - 165 * 0.0245,
  1925 - 165 * 0.0289,
  3569 - 165 * 0.0421
)

X <- solve(A, B)
names(X) <- c("Procedure Cost (P)", "Annual ICD Cost (D)", "Annual No-ICD Cost (C)")
print(X)
