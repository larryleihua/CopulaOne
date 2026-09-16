# Shared validation and numerically stable scalar transforms.
.positive_parameters <- function(...) {
  values <- list(...)
  if (!all(vapply(values, function(x) is.numeric(x) && length(x) == 1L &&
                  is.finite(x) && x > 0, logical(1))))
    stop("Shape parameters must be finite positive scalars.")
}

.positive_quantile <- function(u, cdf, maxit) {
  .probabilities(u)
  if (length(maxit)!=1L || !is.finite(maxit) || maxit<1 || maxit!=floor(maxit))
    stop("maxit must be a positive integer.")
  vapply(u, function(p) {
    if (p==0) return(0)
    if (p==1) return(Inf)
    residual <- function(t) {
      value <- cdf(exp(t))
      if (length(value)!=1L || !is.finite(value) || value<0 || value>1)
        stop("Invalid marginal CDF during quantile inversion.")
      value-p
    }
    lo <- -1; hi <- 1
    while (residual(lo)>0 && lo > -700) lo <- max(2*lo,-700)
    while (residual(hi)<0 && hi < 700) hi <- min(2*hi,700)
    if (residual(lo)>0 || residual(hi)<0)
      stop("Quantile exceeds the supported floating-point range.")
    root <- uniroot(residual,c(lo,hi),tol=1e-10,maxiter=maxit,check.conv=TRUE)$root
    if (abs(residual(root))>max(1e-12,1e-7*min(p,1-p)))
      stop("Quantile inversion did not achieve the required probability accuracy.")
    exp(root)
  },numeric(1),USE.NAMES=FALSE)
}

.probabilities <- function(x, interior = FALSE, name = "probabilities") {
  if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0 | x > 1) ||
      (interior && any(x == 0 | x == 1)))
    stop(name, if (interior) " must be finite and strictly in (0,1)." else
      " must be finite and in [0,1].")
  x
}

.sample_size <- function(n) {
  if (!is.numeric(n) || length(n)!=1L || !is.finite(n) || n<0 ||
      n>.Machine$integer.max || n!=floor(n)) stop("n must be a nonnegative integer.")
  as.integer(n)
}

# Gamma(a) = Gamma(a+1) U^(1/a) for a<1, evaluated on the log scale.
# This avoids zero-draw retry loops and their implicit tail truncation.
.log_gamma_draw <- function(n, shape) {
  .positive_parameters(shape)
  if(shape<1) log(rgamma(n,shape+1)) + log(pmax(runif(n),.Machine$double.xmin))/shape else
    log(rgamma(n,shape))
}

.quadrature_rule <- function(rule) {
  if (!is.list(rule) || !is.numeric(rule$nodes) || !is.numeric(rule$weights) ||
      length(rule$nodes)<2L || length(rule$nodes)!=length(rule$weights) ||
      any(!is.finite(rule$nodes) | rule$nodes<=0 | rule$nodes>=1) ||
      any(!is.finite(rule$weights) | rule$weights<=0) || abs(sum(rule$weights)-1)>1e-6)
    stop("Quadrature requires at least two finite interior nodes and positive weights summing to one.")
  rule
}

.checked_cubature <- function(f, ..., lowerLimit, upperLimit) {
  result <- cubature::adaptIntegrate(f,...,lowerLimit=lowerLimit,upperLimit=upperLimit,
                                     tol=1e-6,absError=1e-9,maxEval=1000000)
  if(length(result$integral)!=1L || !is.finite(result$integral) ||
      length(result$error)!=1L || !is.finite(result$error) ||
      (!is.null(result$returnCode) && result$returnCode!=0L) ||
      result$error>max(1e-9,1e-6*abs(result$integral)))
    stop("Dependence integration did not meet its requested accuracy.")
  result
}

.log1pexp <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
.logexpm1 <- function(x) {
  out <- x + log(-expm1(-x))
  out[x == 0] <- -Inf
  out
}
.logcoshm1 <- function(x) {
  out <- x - log(2) + 2 * log(-expm1(-x))
  out[x < 1e-5] <- 2 * log(x[x < 1e-5]) - log(2) + x[x < 1e-5]^2 / 12
  out
}
.acosh_exp <- function(lx) {
  out <- lx + log(2)
  regular <- lx < 20
  out[regular] <- acosh1p_FRA1(exp(lx[regular]))
  out
}

# F, F', and -F'' evaluated without constructing large powers of s.
.fra_F <- function(ls, p) {
  p <- abs(p)
  A <- .acosh_exp(ls)
  ld <- .5 * (ls + .log1pexp(ls - log(2)) + log(2))
  if (p < 1e-7) {
    value <- 2 * log(A) - log(2)
    first <- log(A) - ld
    ratio <- 1 + 2 / expm1(2 * A) - 1 / A
  } else {
    value <- .logcoshm1(p * A) - 2 * log(p)
    first <- p * A - log(2) + log(-expm1(-2 * p * A)) - log(p) - ld
    ratio <- (1 - p) + 2 / expm1(2 * A) - 2 * p / expm1(2 * p * A)
  }
  second <- first + log(pmax(ratio, 0)) - ld
  small <- ls < log(1e-5)
  if (any(small)) {
    s <- exp(ls[small])
    value[small] <- ls[small] + log1p(-(1-p^2)*s/6 + (1-p^2)*(4-p^2)*s^2/90)
    first[small] <- log1p(-(1-p^2)*s/3 + (1-p^2)*(4-p^2)*s^2/30)
    second[small] <- log((1-p^2)/3) + log1p(-(4-p^2)*s/5)
  }
  if (p == 1) {
    value <- ls; first <- rep(0, length(ls)); second <- rep(-Inf, length(ls))
  }
  list(value = value, first = first, second = second)
}

.fra_Finv <- function(lx, p) {
  p <- abs(p)
  A <- if (p < 1e-7) exp((lx + log(2))/2) else .acosh_exp(lx + 2*log(p))/p
  out <- .logcoshm1(A)
  small <- lx < -20
  out[small] <- lx[small] + log1p((1-p^2)/6 * exp(lx[small]))
  out
}

.fra_J <- function(lt, p) {
  if (p == 1) return(list(value = lt, first = rep(0, length(lt)),
                          second = rep(-Inf, length(lt))))
  lx <- -lt
  ld <- .fra_Finv(0, p)
  fx <- .fra_F(lx, p)
  lT <- .fra_Finv(.log1pexp(fx$value), p)
  fT <- .fra_F(lT, p)
  lS <- lT + log(-expm1(ld - lT))
  lSp <- fx$first - fT$first
  # x S''/S' is a difference of two dimensionless terms.
  xSpp_Sp <- -exp(lx + fx$second - fx$first) +
    exp(lx + fx$first + fT$second - 2*fT$first)
  curvature <- 2 + xSpp_Sp - 2 * exp(lx + lSp - lS)
  small <- lx < log(1e-5)
  if (any(small)) {
    fd <- .fra_F(ld, p)
    s1 <- exp(-fd$first)
    s2 <- -(1-p^2)/3 * s1 + exp(fd$second - 3*fd$first)
    x <- exp(lx[small])
    lS[small] <- lx[small] + log(s1) + log1p(.5*s2/s1*x)
    lSp[small] <- log(s1) + log1p(s2/s1*x)
    # J(t) is asymptotically linear. Its residual curvature is negligible
    # at this order; avoid cancellation of the three O(1) terms above.
    curvature[small] <- 0
  }
  if (any(curvature < -1e-7, na.rm = TRUE))
    stop("FRA1 derivative lost numerical precision.")
  first <- lSp - 2*lt - 2*lS
  list(value = -lS, first = first,
       second = first - lt + log(pmax(curvature, 0)))
}

.fra_Jinv <- function(ly, p) {
  if (p == 1) return(ly)
  ld <- .fra_Finv(0, p)
  lz <- -ly
  lf <- .fra_F(logspace_add_FRA1(ld, lz), p)$value
  delta <- .logexpm1(pmax(lf, 0))
  small <- lz < log(1e-5)
  if (any(small)) {
    fd <- .fra_F(ld, p)
    delta[small] <- fd$first + lz[small] +
      log1p(-.5 * exp(fd$second-fd$first+lz[small]))
  }
  -.fra_Finv(delta, p)
}

.fra_W <- function(lt, theta) {
  if (theta <= 0) return(.fra_J(lt, -theta))
  a <- 1-theta
  h <- .log1pexp(a*lt)
  lk <- .logexpm1(h/a)
  # log1p(exp(a*lt)) can underflow; K(t) ~ t^a/a there.
  tiny <- a*lt < -700
  lk[tiny] <- a*lt[tiny] - log(a)
  kp <- (a-1)*lt + (1/a-1)*h
  kpp <- log1p(-a) + (a-2)*lt + (1/a-2)*h
  j <- .fra_J(lk, 0)
  list(value = j$value, first = j$first+kp,
       second = logspace_add_FRA1(j$second+2*kp, j$first+kpp))
}

.fra_logphi <- function(u, eta, theta) {
  lx <- log(-log(u))
  if (eta > 0) {
    b <- b_FRA1(eta)
    lx <- log(b) + .logexpm1(-log(u)/b)
  }
  ly <- .fra_Finv(lx, abs(eta))
  lj <- .fra_Jinv(ly, if (theta <= 0) -theta else 0)
  if (theta <= 0) return(lj)
  a <- 1-theta
  out <- .logexpm1(a*.log1pexp(lj))/a
  tiny <- lj < -700
  out[tiny] <- (log(a)+lj[tiny])/a
  out
}

.fra_logderiv <- function(lt, eta, theta) {
  w <- .fra_W(lt, theta)
  f <- .fra_F(w$value, abs(eta))
  if (eta <= 0) {
    lg <- -exp(f$value)
    rate <- f$first
    curvature <- logspace_add_FRA1(2*f$first, f$second)
  } else {
    b <- b_FRA1(eta)
    lq <- .log1pexp(f$value-log(b))
    lg <- -b*lq
    rate <- f$first-lq
    curvature <- logspace_add_FRA1(f$second-lq, log1p(1/b)+2*rate)
  }
  list(value = lg, log_negative_first = lg+rate+w$first,
       log_second = lg+logspace_add_FRA1(curvature+2*w$first, rate+w$second))
}
