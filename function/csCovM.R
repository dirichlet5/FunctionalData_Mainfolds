csCovM = function (yList, tList, error = TRUE) 
{
  if (!error) {
    yAr <- simplify2array(yList)
    dims <- dim(yAr)
    n <- dims[3]
    yMat <- matrix(yAr, dims[1] * dims[2], dims[3])
    res <- array(tcrossprod(yMat)/(n - 1), c(dims[1], dims[2], 
                                             dims[1], dims[2]))
    res <- aperm(res, c(2, 4, 1, 3))
  }
  else if (error) {
    stop("Cross-sectional covariance estimate with error is not implemented")
  }
  res
}