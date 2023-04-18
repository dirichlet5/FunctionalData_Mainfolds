RFPCA_mu_trim = function (Ly, Lt, optns = list(), mu_trim) 
{
  optns <- SetOptionsRFPCA(Ly, Lt, optns)
  varnames <- rownames(Ly[[1]])
  for (optn in names(optns)) {
    expr <- parse(text = sprintf("%s <- optns[['%s']]", 
                                 optn, optn))
    eval(expr)
  }
  nRegGrid <- dim(mu_trim)[2] #---------------已修改---------
  Ymat <- do.call(cbind, Ly)
  Tvec <- do.call(c, Lt)
  obsGrid <- sort(unique(Tvec))
  dimAmbient <- nrow(Ymat)
  dimTangent <- calcTanDim(mfd, dimAmbient = dimAmbient)
  ord <- order(Tvec)
  regGrid <- seq(min(obsGrid), max(obsGrid), length.out = nRegGrid)
  a <- ifelse(is.infinite(ToutRange[1]), min(obsGrid), ToutRange[1])
  b <- ifelse(is.infinite(ToutRange[2]), max(obsGrid), ToutRange[2])
  workGrid <- seq(a, b, length.out = nRegGrid)
  m <- length(workGrid)
  if (dataType == "Sparse") {
    if (userBwMu == "GCV") {
      userBwMu <- GCVFrechetMeanCurve(mfd, Ly, Lt, kernel, 
                                      0, "Sparse")$bOpt
      userBwCov <- 2 * userBwMu
    }
    muReg <- frechetMeanCurve(mfd, userBwMu, kernel, xin = Tvec[ord], 
                              yin = Ymat[, ord], xout = regGrid, npoly = npoly)
    muObs <- apply(muReg, 1, function(x) {
      ConvertSupport(as.numeric(regGrid), as.numeric(obsGrid), 
                     mu = x)
    })
    muObs <- apply(muObs, 1, project, mfd = mfd)
  }
  else if (dataType == "Dense") {
    muObs <- t(plyr::daply(data.frame(x = Tvec, t(Ymat)), 
                           "x", function(dat) {
                             y <- t(as.matrix(dat[-1]))
                             c(frechetMean(mfd, y))
                           }, .drop_o = FALSE))
    if (obsGridOnly) {
      muReg <- muObs
    }
    else {
      muReg <- t(apply(muObs, 1, function(x) {
        ConvertSupport(as.numeric(obsGrid), as.numeric(regGrid), 
                       mu = x)
      }))
      muReg <- apply(muReg, 2, project, mfd = mfd)
    }
  }
  if (all.equal(regGrid, workGrid)) {
    muWork <- muReg
  }
  else {
    muWork <- apply(muReg, 1, function(x) {
      ConvertSupport(as.numeric(regGrid), workGrid, mu = x)
    })
    muWork <- apply(muWork, 1, project, mfd = mfd)
  }
  rownames(muObs) <- varnames
  colnames(muObs) <- obsGrid
  rownames(muWork) <- varnames
  colnames(muWork) <- workGrid
  if (optns$meanOnly == TRUE) {
    return(list(muReg = muReg, regGrid = regGrid, muWork = muWork, 
                workGrid = workGrid))
  }
  if (meanOnly) {
    res <- list(muReg = muReg, muWork = muWork, muObs = muObs, 
                regGrid = regGrid, workGrid = workGrid, obsGrid = obsGrid, 
                userBwMu = userBwMu, userBwCov = userBwCov, mfd = mfd, 
                optns = optns)
    class(res) <- "RFPCA"
    return(res)
  }
  yListLog <- lapply(seq_along(Ly), function(i) {
    tt <- Lt[[i]]
    mu0 <- muObs[, match(tt, obsGrid), drop = FALSE]
    yy <- rieLog(mfd, mu_trim, Ly[[i]])
    yy
  })
  if (dataType == "Sparse") {
    covWork <- smoothCovM2(yListLog, Lt, matrix(0, nrow(muObs), 
                                                ncol(muObs)), workGrid, userBwCov, kernel, error = error)
  }
  else if (dataType == "Dense") {
    covObs <- csCovM(yListLog, Lt, error = FALSE)
    if (obsGridOnly) {
      covWork <- covObs
    }
    else {
      covWork <- array(apply(covObs, c(3, 4), function(mat) {
        ConvertSupport(as.numeric(obsGrid), as.numeric(workGrid), 
                       Cov = mat, isCrossCov = TRUE)
      }), c(length(workGrid), length(workGrid), dimTangent, 
            dimTangent))
    }
  }
  covWork <- projectCov(mfd, covWork, muWork)
  tmp <- lapply(seq_along(Lt), function(i) {
    tt <- Lt[[i]]
    yy <- yListLog[[i]]
    ind <- tt >= ToutRange[1] & tt <= ToutRange[2]
    list(t = tt[ind], y = yy[, ind, drop = FALSE])
  })
  LtTrunc <- sapply(tmp, `[[`, "t", simplify = FALSE)
  LyLogTrunc <- sapply(tmp, `[[`, "y", simplify = FALSE)
  names(LtTrunc) <- names(LyLogTrunc) <- names(Lt)
  subInd <- sapply(LtTrunc, length) > 0
  if (!all(subInd)) {
    warning("Some subjects have no observations!")
    LtTrunc <- LtTrunc[subInd]
    LyLogTrunc <- LyLogTrunc[subInd]
  }
  TTruncInd <- obsGrid >= ToutRange[1] & obsGrid <= ToutRange[2]
  obsGridTrunc <- obsGrid[TTruncInd]
  mObsTrunc <- length(obsGridTrunc)
  muObsTrunc <- muObs[, TTruncInd, drop = FALSE]
  if (!is.null(userSigma2)) {
    sigma2 <- userSigma2
  }
  else if (error) {
    sigma2 <- EstSigma2(mfd, LyLogTrunc, LtTrunc, matrix(0, 
                                                         nrow(LyLogTrunc[[1]]), ncol(muObsTrunc)), covWork, 
                        workGrid, userBwCov, kernel, smooth = FALSE)
  }
  else {
    sigma2 <- 0
  }
  eig <- EigenAnalysis(covWork, workGrid, maxK, FVEthreshold, 
                       fastEig = fastEig, verbose = verbose)
  functionalObs <- inherits(mfd, "L2")
  if (functionalObs) {
    grids <- 1/(dimTangent - 1)
    eig[["lambda"]] <- eig[["lambda"]] * grids
    eig[["phi"]] <- eig[["phi"]]/sqrt(grids)
  }
  lam <- eig[["lambda"]]
  K <- length(lam)
  phi <- array(eig[["phi"]], c(m, dimTangent, K))
  covFitted <- eig[["fittedCov"]]
  if (methodMuCovEst == "smooth") {
    covObs <- apply(eig[["fittedCov"]], c(3, 4), function(mat) {
      ConvertSupport(as.numeric(workGrid), as.numeric(obsGridTrunc), 
                     Cov = mat, isCrossCov = TRUE)
    })
    covObs <- array(covObs, c(mObsTrunc, mObsTrunc, dimTangent, 
                              dimTangent))
    covObs <- projectCov(mfd, covObs, muObsTrunc)
  }
  if (obsGridOnly) {
    phiObsTrunc <- phi
  }
  else {
    phiObsTrunc <- apply(phi, c(2, 3), function(phi) {
      ConvertSupport(as.numeric(workGrid), as.numeric(obsGridTrunc), 
                     mu = phi)
    })
  }
  phiObsTrunc <- projectPhi(mfd, phiObsTrunc, muObsTrunc)
  if (methodXi == "CE") {
    CE <- CEScores(LyLogTrunc, LtTrunc, list(), matrix(0, 
                                                       dimTangent, mObsTrunc), obsGridTrunc, covObs, lam, 
                   phiObsTrunc, sigma2)
    xi <- t(do.call(cbind, CE["xiEst", ]))
    if (functionalObs) {
      stop("Not implemented yet")
    }
  }
  else if (methodXi == "IN") {
    ylogTrunc <- aperm(simplify2array(LyLogTrunc), c(3, 
                                                     2, 1))
    dims <- dim(ylogTrunc)
    n <- dims[1]
    yMat2 <- matrix(ylogTrunc, n, dims[2] * dims[3])
    phiMat <- matrix(phiObsTrunc, ncol = dim(phiObsTrunc)[3])
    if (length(obsGridTrunc) == 1) {
      gridt <- 1
    }
    else {
      gridt <- mean(diff(obsGridTrunc))
    }
    xi <- yMat2 %*% phiMat * gridt
    if (functionalObs) {
      xi <- xi * grids
    }
  }
  rownames(xi) <- names(LtTrunc)
  colnames(xi) <- paste0("xi", seq_len(ncol(xi)))
  res <- list(muReg = muReg, muWork = muWork, muObs = muObs, 
              muObsTrunc = muObsTrunc, cov = covWork, covFitted = covFitted, 
              phi = phi, covObs = covObs, phiObsTrunc = phiObsTrunc, 
              lam = lam, xi = xi, sigma2 = sigma2, regGrid = regGrid, 
              workGrid = workGrid, obsGrid = obsGrid, obsGridTrunc = obsGridTrunc, 
              K = K, userBwMu = userBwMu, userBwCov = userBwCov, mfd = mfd, 
              optns = optns)
  class(res) <- "RFPCA"
  res
}
