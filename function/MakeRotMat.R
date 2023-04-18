MakeRotMat = function(p1, p2, tol=1e-10) {
  d <- length(p1)
  if (length(p1) != length(p2)) {
    stop('dimensions of p1 and p2 must agree')
  }
  if (sum(abs(p1 - Normalize(p1))) > tol || 
      sum(abs(p2 - Normalize(p2))) > tol) {
    stop('p1 and p2 needs to be have norm 1')
  }
  
  if (sum(abs(p1 + p2)) < tol) {
    warning('Rotation between podal points is arbitrary')
  }
  
  Sphere <- structure(1, class='Sphere')
  psi <- as.numeric(distance(Sphere, matrix(p1), matrix(p2)))
  p1t <- Normalize(Proj(p1, p2, rej=TRUE))
  R <- diag(1, nrow=d) + 
    sin(psi) * (tcrossprod(p2, p1t) - tcrossprod(p1t, p2)) + 
    (cos(psi) - 1) * (tcrossprod(p2, p2) + tcrossprod(p1t, p1t))
  
  R
}
