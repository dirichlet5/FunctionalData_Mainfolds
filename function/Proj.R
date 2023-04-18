Proj = function(p1, p2, rej=FALSE, tol=1e-10) {
  p1 <- as.numeric(p1)
  p2 <- Normalize(as.numeric(p2), tol)
  
  if (sum(p2^2) == 0) {
    stop('p2 cannot be 0')
  }
  
  proj <- c(crossprod(p1, p2)) * p2
  res <- if (!rej) {
    proj
  } else { # Return the rejection
    p1 - proj
  }
  
  as.numeric(res)
}
