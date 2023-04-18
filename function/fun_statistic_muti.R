fun_statistic_muti <- function(n, m, y_sample, mu_trim, yList) {
  Resfpca <- RFPCA_mu_trim(y_sample$Ly, y_sample$Lt, list(maxK = 5, mfd = mfd), mu_trim)
  
  xi <- matrix(data = 0, nrow = n, ncol = maxK)
  yList_tange = list()
  for (i in 1:n) {
    yList_tange[[i]] <- rieLog.Sphere(mfd, p = mu_trim, yList[[i]])
  }
  for (i in 1:n) {
    for (j in 1:3) {
      temp <- yList_tange[[i]][j, ] %*% Resfpca$phi[, j, ]
      xi[i,] = xi[i,] + temp
    }
    
  }
  T_stat <- (xi^2) %*% (1 / Resfpca$lam)
  
  
  list(T_stat = T_stat)
}
