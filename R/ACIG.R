# ACIG.R

logK <- function(x, nu)
{
  N <- length(x)
  rst <- rep(0,N)
  out <- .C("lnK", as.double(nu), as.double(x), as.integer(N), rst=as.double(rst))
  out$rst
}

logG <- function(alp)
{
  N <- length(alp)
  rst <- rep(0,N)
  out <- .C("lnG", as.double(alp), as.integer(N), rst=as.double(rst))
  out$rst
}

logden_Gumbel <- function(F1, F2, alp, lf1, lf2)
{
  N <- length(F1)
  rst <- rep(0,N)
  out <- .C("logden_Gumbel", as.double(F1), as.double(F2), as.double(alp), as.double(lf1), as.double(lf2), as.integer(N), rst=as.double(rst))
  out$rst
}


## optimized version
igltinv <- function(t,alp)
{
  N <- length(t)
  rst <- rep(0,N)
  out <- .C("invpsi", as.double(t), as.double(alp), as.integer(N), rst <- as.double(rst))
  out$rst
}

logden_ACIG <- function(F1, F2, alp, lf1, lf2)
{
  N <- length(F1)
  s <- s1 <- s2 <- logden <- rep(0, N)
  s1 <- igltinv(F1, alp)
  s2 <- igltinv(F2, alp)
  s <- s1 + s2
  logden <- logK(2*sqrt(s), alp-2) - logK(2*sqrt(s1), alp-1) - logK(2*sqrt(s2), alp-1) + ((alp -2)/2)*log(s) - ((alp-1)/2)*(log(s1)+log(s2)) + lf1 + lf2 - log(2) + logG(alp)
  logden
}

iglt <- function(s,alp)
{
  N <- length(s)
  rst <- rep(0,N)
  out <- .C("psi", as.double(s), as.double(alp), as.integer(N), rst <- as.double(rst))
  out$rst
}  

rACIG <- function(nsim,alp)
{
  if(alp<0) stop("alp>=0")
  uvec <- rep(0,nsim)
  vvec <- rep(0,nsim)
  # sample inverse Gamma V
  V <- rinvgamma(nsim, shape <- alp, scale <- 1)
  X1 <- runif(nsim)
  X2 <- runif(nsim)
  U1 <- iglt(-log(X1)/V, rep(alp,nsim))
  U2 <- iglt(-log(X2)/V, rep(alp,nsim))
  cbind(U1, U2)
}
