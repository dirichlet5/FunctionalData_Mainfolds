outlier_detect <- function(samp, outlier_pos) {
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

  distance_mu <- c()
  for (i in 1:n) {
    d_M_p1p2 <- distance(mfd, samp$X[i, 1:3, ], muMt)
    distance_mu_temp <- t(d_M_p1p2) %*% d_M_p1p2 # 距离的模长
    distance_mu <- c(distance_mu, distance_mu_temp)
  }
  # 选出距离最近的80%的样本
  mu_trim_pos <- order(distance_mu, decreasing = FALSE)[1:(n * 0.8)]
  # 计算筛选后样本的均值mu_trim
  spSamp_trim <- SparsifyM(samp$X[mu_trim_pos, , ], samp$T, sparsity)
  yList_trim <- spSamp_trim$Ly
  tList_trim <- spSamp_trim$Lt
  mu_trim <- RFPCA(yList_trim, tList_trim, list(maxK = K, mfd = mfd))$muObs

  source("./function/RFPCA_mu_trim.R")
  source("./function/SetOptionsRFPCA.R")
  source("./function/csCovM.R")
  source("./function/projectCov.R")
  source("./function/EigenAnalysis.R")
  source("./function/projectPhi.R")

  Resfpca <- RFPCA_mu_trim(yList, tList, list(maxK = maxK, mfd = mfd), mu_trim)

  xi <- matrix(data = 0, nrow = n, ncol = maxK)
  # 将yList映射到切空间上，记为yList_tange
  yList_tange <- list()
  for (i in 1:n) {
    yList_tange[[i]] <- rieLog.Sphere(mfd, p = mu_trim, yList[[i]])
  }
  for (i in 1:n) {
    for (j in 1:3) {
      #  (yList需要映射到切空间上，再和phi做内积)
      temp <- yList_tange[[i]][j, ] %*% Resfpca$phi[, j, ]
      xi[i, ] <- xi[i, ] + temp
    }
  }

  T_stat <- (xi^2) %*% (1 / Resfpca$lam)

  source("./function/fun_statistic_muti.R")

  k_sub <- 0.75

  iter <- 2
  T_matrix <- matrix(0, iter, n)
  for (j in 1:iter) {
    ind <- sample(c(1:n), floor(n / 2), replace = FALSE)
    y_sample <- list()

    for (k in 1:length(ind)) {
      y_sample$Lt[[k]] <- spSamp$Lt[[ind[k]]]
      y_sample$Ly[[k]] <- spSamp$Ly[[ind[k]]]
    }

    T_matrix[j, ] <- fun_statistic_muti(n, m, y_sample, mu_trim, yList)$T_stat
  }

  T_max <- apply(T_matrix, 2, max)
  T_min <- apply(T_matrix, 2, min)
  clean_max <- order(T_max, decreasing = FALSE)[1:floor(n * k_sub)]
  clean_min <- order(T_min, decreasing = FALSE)[1:floor(n * k_sub)]
  cleanset <- intersect(clean_max, clean_min)

  # 记录clean_set
  y_clean <- list()
  for (i in 1:length(cleanset)) {
    y_clean$Lt[[i]] <- spSamp$Lt[[cleanset[i]]]
    y_clean$Ly[[i]] <- spSamp$Ly[[cleanset[i]]]
  }

  # 计算Q_i
  # Q_i(根据T_i计算的值，也是一维的)
  Q_T <- rep(0, n)

  T.test <- fun_statistic_muti(n, m, y_sample, mu_trim, yList)$T_stat
  # T.test <- fun_statistic_muti(n, m, spSamp, mu_trim, yList)$T_stat

  th_median <- median(T.test) / qchisq(0.5, maxK)


  T.test <- T.test / th_median
  # print(T.test)


  for (i in 1:n) {
    Q_T[i] <- qnorm((pchisq(T.test[i], maxK) + 1) / 2) # 计算Q_i
  }

  len.tt <- n
  U_thr <- sqrt(2 * log(n))
  t_it <- seq(0, U_thr, length = len.tt)
  t_fdr <- rep(0, len.tt)
  for (it in 1:len.tt) {
    t_fdr[it] <- 2 * n * (1 - pnorm(t_it[it])) / max(length(which(abs(Q_T) >= t_it[it])), 1)
  }

  # 计算t0，# 添加一个偏移量
  alpha0 <- 0.05
  if (length(which(t_fdr <= alpha0)) > 0) {
    t0 <- t_it[min(which((t_fdr - alpha0) <= 0))]
  } else {
    t0 <- sqrt(2 * log(n))
  }

  # paste("t0:", t0)

  # 检验大于t0的为outlier
  out_set <- which(Q_T > t0)

  source("./function/evaluation.R")

  tpr <- evaluation(n = n, results = out_set, true = outlier_pos)$tp
  fpr <- evaluation(n = n, results = out_set, true = outlier_pos)$fp
  list(tpr = tpr, fpr = fpr)
}
