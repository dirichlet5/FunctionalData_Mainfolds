MakeSphericalProcess = function(
  n, mu, pts=seq(0, 1, length.out=ncol(mu)), K=2, lambda=rep(1, K), sigma2=0, 
  xiFun = rnorm, epsFun = rnorm, epsBase = c('mu', 'X'), 
  basisType=c('cos', 'sin', 'fourier', 'legendre01')) {
  
  epsBase <- match.arg(epsBase)
  basisType <- match.arg(basisType)
  p <- nrow(mu)
  m <- length(pts)
  
  # A function that generates iid centered and scaled random variables for scores
  xi <- matrix(xiFun(n * K), n, K) %*% diag(sqrt(lambda), nrow = K)
  phi <- matrix(MakephiS(K, mu, pts, basisType), ncol=K)
  xiphi <- array(xi %*% t(phi), c(n, m, p))
  
  mfd <- structure(1, class='Sphere')
  X <- sapply(seq_len(m), function(tt) {
    mu0 <- mu[, tt]
    V <- matrix(xiphi[, tt, ], dim(xiphi)[1], dim(xiphi)[3])
    t(rieExp(mfd, mu0, t(V)))
  }, simplify='array')
  dimnames(X) <- list(i = seq_len(n), j = seq_len(p), t = pts)
  
  if (sigma2 > 1e-15) {
    if (epsBase == 'mu') {
      epsArray <- vapply(seq_len(m), function(tt) {
        rot <- MakeRotMat(c(rep(0, p - 1), 1), mu[, tt])
        eps <- cbind(matrix(epsFun(n * (p - 1)), n, p - 1), 0) %*% t(rot) * 
          sqrt(sigma2)
        eps
      }, matrix(0, n, p))
      xiphiN <- xiphi + aperm(epsArray, c(1, 3, 2))
      XNoisy <- sapply(seq_len(m), function(tt) {
        mu0 <- mu[, tt]
        V <- matrix(xiphiN[, tt, ], dim(xiphiN)[1], dim(xiphiN)[3])
        t(rieExp(mfd, mu0, t(V)))
      }, simplify='array')
    } else if (epsBase == 'X') {
      XNoisy <- apply(X, c(1, 3), function(x) {
        rot <- MakeRotMat(c(rep(0, p - 1), 1), x)
        eps <- rot %*% matrix(c(epsFun(p - 1), 0) * sqrt(sigma2))
        t(rieExp(mfd, x, eps))
      })
      XNoisy <- aperm(XNoisy, c(2, 1, 3))
    }
  } else {
    XNoisy <- X
  }
  
  
  res <- list(XNoisy = XNoisy, X = X, T = pts, xi = xi, phi=phi)
  res
}
