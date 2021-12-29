n <- 1000
mat_size <- 1000
mats <- list()
start_time <- NULL
end_time <- NULL
for (i in 1:n)
{
  cat("last iteration took", as.numeric(end_time - start_time), "secs\n")
  start_time <- Sys.time()
  cat("iteration", i, "running")
  a_mat <- matrix(rnorm(mat_size^2), nrow = mat_size, ncol=mat_size)
  b_mat <- matrix(rnorm(mat_size^2), nrow = mat_size, ncol=mat_size)
  mats[[i]] <- solve(a_mat %*% a_mat)
  cat('\014')
  end_time <- Sys.time()
}
