ord_detect <- function(samp, outlier_pos) {
  n <- dim(samp$X)[1]
  m <- dim(samp$X)[3]
  sparsity <- m
  spSamp <- SparsifyM(samp$X, samp$T, sparsity)
  yList <- spSamp$Ly
  tList <- spSamp$Lt

  maxK <- 5
  resSp <- RFPCA(yList, tList, list(maxK = maxK, mfd = mfd))
  # resSp$lam # 取得的maxK要大于0.001
  muMt <- resSp$muObs

  # 计算distance
  distance_mu <- c()
  for (i in 1:n) {
    d_M_p1p2 <- distance(mfd, samp$X[i, 1:3, ], muMt)
    distance_mu_temp <- t(d_M_p1p2) %*% d_M_p1p2 # 距离的模长
    distance_mu <- c(distance_mu, distance_mu_temp)
  }

  # 根据distance选择一个阈值进行拒绝
  # 选出距离最远的20个为异常值
  mu_trim_pos <- order(distance_mu, decreasing = TRUE)[1:(n * 0.05)] # 阈值0.05

  source("./function/evaluation.R")

  tpr <- evaluation(n = n, results = mu_trim_pos, true = outlier_pos)$tp
  fpr <- evaluation(n = n, results = mu_trim_pos, true = outlier_pos)$fp
  list(tpr = tpr, fpr = fpr)
}
