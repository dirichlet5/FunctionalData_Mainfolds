mfd <- structure(1, class = "Sphere")

rieExp.Sphere <- function(mfd, p, V) {
  tol <- 1e-10

  V <- as.matrix(V)
  p <- as.matrix(p)
  m <- dim(p)[2]
  n <- dim(V)[2]
  d <- dim(V)[1]

  if (length(p) == 0 || length(V) == 0) {
    return(matrix(0, nrow(p), 0))
  }

  if (m == 1 && n != 1) {
    p <- matrix(p, d, n)
  } else if (m > 1 && n == 1) {
    V <- matrix(V, d, m)
  }
  mn <- max(m, n)

  # p needs to be on a unit sphere
  # stopifnot(abs(sum(p^2) - 1) <= tol)

  stopifnot(nrow(V) == nrow(p))

  # each col of V needs to be orthogonal to p
  stopifnot(all(abs(colSums(p * V)) <= tol))

  res <- vapply(seq_len(mn), function(i) {
    Exp1(V[, i], p[, i], tol = tol)
  }, rep(0, d))
  if (!is.matrix(res)) {
    res <- matrix(res, nrow = 1)
  }
  res
}

Exp1 <- function(v, mu, tol = 1e-10) {
  vNorm <- as.numeric(sqrt(crossprod(v)))
  if (!is.na(vNorm) && vNorm <= tol) {
    mu
  } else {
    cos(vNorm) * mu + sin(vNorm) * v / vNorm
  }
}

rieLog.Sphere <- function(mfd, p, X, tol = 1e-10) {
  X <- as.matrix(X)
  p <- as.matrix(p)
  m <- dim(p)[2]
  n <- dim(X)[2]
  d <- dim(X)[1]

  if (length(p) == 0 || length(X) == 0) {
    return(matrix(0, calcTanDim.Sphere(dimAmbient = d), 0))
  }

  if (m == 1 && n != 1) {
    p <- matrix(p, d, n)
  } else if (m > 1 && n == 1) {
    X <- matrix(X, d, m)
  }
  mn <- max(m, n)

  # # p needs to be on a unit sphere
  # stopifnot(all(abs(colSums(p^2) - 1) <= tol))

  # # each column of X needs to be on a unit sphere
  # stopifnot(all(abs(apply(X, 2, function(x) sum(x^2)) - 1) <= tol))

  # Z <- vapply(seq_len(n), function(i) {
  # Log1(X[, i], p[, i], tol)
  # }, rep(0, d))
  Z <- Log2(X, p, tol)

  return(Z)
}

Log2 <- function(X, Mu, tol = 1e-10) {
  dimAmbient <- nrow(X)
  N <- ncol(X)
  cprod <- colSums(X * Mu)
  U <- X - matrix(cprod, dimAmbient, N, byrow = TRUE) * Mu
  uNorm <- sqrt(colSums(U^2))
  res <- matrix(0, dimAmbient, N)
  ind <- uNorm > tol
  distS <- acos(ifelse(cprod > 1, 1, ifelse(cprod < -1, -1, cprod)))
  res[, ind] <- (U * matrix(distS / uNorm, dimAmbient, N, byrow = TRUE))[, ind, drop = FALSE]
  res
}

muList <- list(
  function(x) x * 2,
  function(x) sin(x * 1 * pi) * pi / 2 * 0.6,
  function(x) rep(0, length(x))
)

n <- 500 # 样本量
# outlier_n <- 20 # 异常值的数量
m <- 51
pts <- seq(0, 1, length.out = m)
mu <- Makemu(mfd, muList, c(rep(0, 2), 1), pts) # 生成均值mu
# matplot(pts, t(mu), type = "l", lty = 1) # mu的坐标xyz

K <- 20
maxK <- 5
sigma2 <- 0

data_type2 <- function(para_k) {
  library("RFPCA")
  library("MASS")
  source("./function/MakephiS.R")
  source("./function/MakeRotMat.R")
  source("./function/MakeSphericalProcess.R")
  source("./function/MakeSphericalOutlier.R")
  source("./function/MakeSphericalOutlier_2.R")
  source("./function/Proj.R")


  # 生成样本数据
  samp <- MakeSphericalProcess(n, mu, pts,
    K = K,
    lambda = 0.07^(seq_len(K) / 2),
    basisType = "legendre01",
    sigma2 = sigma2
  )

  # 异常添加方式2（更换均值函数mu）：
  muList_out <- list(
    function(x) x * 2 * para_k,
    function(x) sin(x * 1 * pi) * pi / 2 * 0.6,
    function(x) rep(0, length(x))
  )
  mu_out <- Makemu(mfd, muList_out, c(rep(0, 2), 1), pts) # 生成均值mu

  outlier_2 <- MakeSphericalProcess(outlier_n, mu_out, pts,
    K = K,
    lambda = 0.07^(seq_len(K) / 2),
    basisType = "legendre01",
    sigma2 = sigma2
  )$X

  outlier_pos <- 1:outlier_n
  samp$X[outlier_pos, , ] <- outlier_2
  return(samp)
}
