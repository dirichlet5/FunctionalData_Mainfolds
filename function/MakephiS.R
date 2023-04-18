MakephiS = function(K, mu, pts=seq(0, 1, length.out=ncol(mu)), 
         type=c('cos', 'sin', 'fourier', 'legendre01')) {
  
  type <- match.arg(type)
  if (!is.matrix(mu)) {
    stop('mu has to be a matrix')
  }
  p <- nrow(mu)
  if (p <= 1) {
    stop('mu must have more than 1 rows')
  }
  m <- length(pts)
  if (ncol(mu) != m) {
    stop('mu must have the same number of columns as the length of pts')
  }
  
  if (min(pts) < 0 || max(pts) > 1) {
    stop('The range of pts should be within [0, 1]')
  }
  
  ptsSqueeze <- do.call(c, lapply(seq_len(p - 1), function(i) {
    pts + (max(pts) + mean(diff(pts))) * (i - 1)
  }))
  ptsSqueeze <- ptsSqueeze / max(ptsSqueeze)
  
  phi1 <- rbind(fdapace:::CreateBasis(K, ptsSqueeze, type),
                matrix(0, m, K)) / sqrt(p - 1) # at the north pole c(0, ..., 0, 1)
  phi1 <- array(phi1, c(m, p, K))
  dimnames(phi1) <- list(t=pts, j=seq_len(p), k=seq_len(K))
  
  # Rotate to the locations of mu
  phi <- sapply(seq_along(pts), function(i) {
    R <- MakeRotMat(c(rep(0, p - 1), 1), mu[, i])
    R %*% phi1[i, , , drop=TRUE]
  }, simplify=FALSE)
  phi <- do.call(abind::abind, list(phi, along=0))
  dimnames(phi) <- dimnames(phi1)
  
  phi
}
