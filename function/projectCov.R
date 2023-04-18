projectCov = function (mfd = structure(1, class = "Sphere"), covR, mu) 
{
  stopifnot(dim(covR)[1] == dim(covR)[2] && dim(covR)[1] == 
              ncol(mu))
  m <- dim(covR)[1]
  projM <- lapply(seq_len(ncol(mu)), function(i) projectTangent(mfd, 
                                                                mu[, i], projMatOnly = TRUE))
  for (j1 in seq_len(m)) {
    for (j2 in seq_len(m)) {
      covR[j1, j2, , ] <- as.matrix(projM[[j1]] %*% covR[j1, 
                                                         j2, , ] %*% Matrix::t(projM[[j2]]))
    }
  }
  covR
}