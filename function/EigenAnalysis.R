EigenAnalysis = function (smoothCov, regGrid, K, FVE = 1, muWork = NULL, tol = 1e-14, 
          verbose = TRUE, fastEig = FALSE) 
{
  if (missing(K)) {
    K <- Inf
  }
  else {
    stopifnot(K >= 1)
  }
  if (length(regGrid) == 1) {
    gridSize <- 1
  }
  else {
    gridSize <- regGrid[2] - regGrid[1]
  }
  nT <- dim(smoothCov)[1]
  p <- dim(smoothCov)[3]
  mat <- matrix(aperm(smoothCov, c(1, 3, 2, 4)), nT * p, nT * 
                  p)
  mat <- (mat + t(mat))/2
  if (fastEig) {
    eig <- tryCatch({
      suppressWarnings(RSpectra::eigs_sym(mat, min(K, 
                                                   nrow(mat)), which = "LA"))
    }, error = function(e) {
      eigen(mat)
    })
  }
  else {
    eig <- eigen(mat)
  }
  positiveInd <- eig[["values"]]/(tol + eig[["values"]][1]) > 
    tol
  if (sum(positiveInd) == 0) {
    stop("All eigenvalues are negative. The covariance estimate is incorrect.")
  }
  d <- eig[["values"]][positiveInd]
  eigenV <- eig[["vectors"]][, positiveInd, drop = FALSE]
  if (K > length(d)) {
    if (verbose && !is.infinite(K)) {
      warning(sprintf("Returning only %d < K = %d components with positive eigenvalues.\n", 
                      length(d), K))
    }
    K <- length(d)
  }
  d <- d[seq_len(K)]
  eigenV <- eigenV[, seq_len(K), drop = FALSE]
  cumFVE <- cumsum(d)/sum(d)
  FVEK <- min(which(cumFVE >= FVE))
  d <- d[seq_len(FVEK)]
  eigenV <- eigenV[, seq_len(FVEK), drop = FALSE]
  if (is.null(muWork)) {
    muWork = seq_len(nrow(eigenV))
  }
  phi <- apply(eigenV, 2, function(x) {
    x <- x/sqrt(gridSize)
    if (0 <= crossprod(x, muWork)) {
      x
    }
    else {
      -x
    }
  })
  lambda <- gridSize * d
  fittedCov <- phi %*% diag(x = lambda, nrow = length(lambda)) %*% 
    t(phi)
  fittedCov <- aperm(array(fittedCov, c(nT, p, nT, p)), c(1, 
                                                          3, 2, 4))
  return(list(lambda = lambda, phi = phi, cumFVE = cumFVE, 
              fittedCov = fittedCov))
}
