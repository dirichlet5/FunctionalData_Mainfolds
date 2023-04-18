projectPhi = function (mfd = structure(1, class = "Sphere"), phi, mu) 
{
  dimTangent <- dim(phi)[2]
  dimAmbient <- nrow(mu)
  m <- dim(phi)[1]
  stopifnot(dimTangent == calcTanDim(mfd, dimAmbient = dimAmbient) && 
              dim(phi)[1] == ncol(mu))
  projM <- lapply(seq_len(ncol(mu)), function(i) projectTangent(mfd, 
                                                                mu[, i], projMatOnly = TRUE))
  for (j in seq_len(m)) {
    phi[j, , ] <- as.matrix(projM[[j]] %*% phi[j, , ])
  }
  phi
}
