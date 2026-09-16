#' Joint cdf of the upper extreme value copula of the GGEE copula
#'
#' Joint cdf of the upper extreme value copula of the GGEE copula
#' @param u Probability or first copula argument.
#' @param v Second copula argument.
#' @param b Dependence parameter; 0 < b < 1 uses the extreme-value formula, otherwise the CDF is u*v.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' pUEV_GGEE_COP(0.3, 0.4, 1.2)
pUEV_GGEE_COP <- function(u, v, b)
{
  if (length(b) != 1L || is.na(b))
    stop("b must be a non-missing scalar.")
  if (b <= 0 || b >= 1) return(u * v)
  if (u == 0 || v == 0) return(0)
  if (u == 1) return(v)
  if (v == 1) return(u)
  if (u == v)
  {
    out <- u^(1 + b)
  } else
  {
    w1 <- -log(u)
    w2 <- -log(v)
    out <- exp(-(w1^(1 + 1/b) - w2^(1 + 1/b))/(w1^(1/b) - w2^(1/b)))
  }
  return(out)
}

pow <- function(x, y){ x^y }
#' Density function of the upper extreme value copula for the GGEE copula
#'
#' Density function of the upper extreme value copula for the GGEE copula
#' @param u Probability or first copula argument.
#' @param v Second copula argument.
#' @param b Dependence parameter; 0 < b < 1 uses the extreme-value formula, otherwise the density is 1.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' dUEV_GGEE_COP(0.3, 0.4, 1.2)
dUEV_GGEE_COP <- function(u, v, b)
{
  if (length(b) != 1L || is.na(b))
    stop("b must be a non-missing scalar.")
  if (b <= 0 || b >= 1) return(1)
  if (u == v)
  {
    out <- (1/12) * (3 * v^b * b^3 * log(v) + 6 * v^b * b^2 * log(v) + 2 * 
      v^b * b^2 + 3 * b * log(v) * v^b - 2 * v^b)/(v * log(v) * b)
  } else
  {
    t1 <- log(u)
    t2 <- 1/b
    t3 <- pow(-t1, t2)
    t4 <- log(v)
    t5 <- pow(-t4, t2)
    t6 <- t3 - t5
    t7 <- 1/t6
    t9 <- pow(u, t7 * t5)
    t11 <- pow(v, -t7 * t3)
    t13 <- 2 * t2
    t14 <- pow(-t1, t13)
    t16 <- (-2 + b) * t2
    t17 <- pow(-t4, -t16)
    t21 <- 3 * t2
    t22 <- pow(-t4, t21)
    t23 <- t3 * t22
    t25 <- pow(-t1, t21)
    t26 <- t25 * t5
    t28 <- pow(-t4, t13)
    t29 <- t14 * t28
    t32 <- (-3 + b) * t2
    t33 <- pow(-t4, -t32)
    t34 <- t3 * t33
    t36 <- pow(-t1, -t32)
    t37 <- t5 * t36
    t39 <- pow(-t1, -t16)
    t44 <- (-1 + b) * t2
    t45 <- pow(-t4, -t44)
    t46 <- t25 * t45
    t47 <- pow(-t1, -t44)
    t48 <- t22 * t47
    t49 <- -2 * t14 * t17 * b + t23 * b + t26 * b - 2 * t29 + t34 * b + t37 * 
      b - 2 * t28 * t39 * b - t46 - t48 + t34 + t37
    t51 <- (2 + b) * t2
    t52 <- pow(-t1, t51)
    t53 <- t52 * t17
    t54 <- pow(-t4, t51)
    t55 <- t39 * t54
    t58 <- (3 + b) * t2
    t59 <- pow(-t4, t58)
    t66 <- b * b
    t68 <- pow(-t1, t58)
    t75 <- t53 + t55 + t55 * b - t47 * t59 * b - 2 * t29 * b + t46 * b + t48 * 
      b + t26 * t66 - t68 * t45 * b + t53 * b - 2 * t29 * t66 + t23 * t66
    t79 <- t6 * t6
    t80 <- t79 * t79
    out <- -t9 * t11 * (t49 + t75)/t66/t80
  }
  return(out)
}

# Shared stable tail-dependence function L(x,y) = max(x,y) H(min/max).
# Combining the old A1/A2 coefficients removes their singularity at
# a - b + be = 0. The A2 coefficient uses b^2, not a^2.
.UEV_PPPP_bundle <- function(w1, w2, be, a, b) {
  wp <- max(w1, w2)
  r <- min(w1, w2) / wp
  k <- b / be
  ell <- 1 + a / be
  P <- be * (a + be) * (a + 2 * b - be) / ((a + b)^2 * (2 * b - be))
  Q <- (a + be) * (b - be) / (a + b)^2
  R <- be * (b - be) * (2 * a + b + be) / ((a + b)^2 * (2 * a + be))

  # D(r) = (r^ell - r^k)/(ell-k).
  # Its exact limit at ell=k is r^k log(r). Factor out the smaller
  # power so expm1 has a nonpositive argument even far into the tails.
  m <- min(k, ell)
  delta <- abs(ell - k)
  log_r <- log(r)
  E <- if (delta == 0) log_r else expm1(delta * log_r) / delta
  z <- exp(delta * log_r)
  rk <- r^k
  rl <- r^ell
  rm <- r^m
  rm1 <- r^(m - 1)
  H <- 1 + P * rk + R * rl - Q * rm * E
  Hp <- P * k * r^(k - 1) + R * ell * r^(ell - 1) -
    Q * rm1 * (m * E + z)
  # Store r H'' rather than H'' to avoid divergent intermediate powers.
  rHpp <- P * k * (k - 1) * r^(k - 1) +
    R * ell * (ell - 1) * r^(ell - 1) -
    Q * rm1 * (m * (m - 1) * E + (2 * m + delta - 1) * z)
  list(value = wp * H, first_min = Hp, first_max = H - r * Hp,
       negative_mixed = rHpp / wp)
}

.check_UEV_PPPP <- function(u, v, be, a, b) {
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!all(vapply(list(be, a, b), scalar, logical(1))) ||
      any(c(be, a, b) <= 0))
    stop("be, a, and b must be finite positive scalars.")
  if (!scalar(u) || !scalar(v) || any(c(u, v) < 0 | c(u, v) > 1))
    stop("Copula arguments must be finite scalar probabilities in [0,1].")
}

A_UEV_PPPP_COP <- function(w1, w2, be, a, b) {
  if (b <= be) return(w1 + w2)
  if (w1 == 0 || w2 == 0) return(max(w1, w2))
  .UEV_PPPP_bundle(w1, w2, be, a, b)$value
}
#' Joint cdf of the upper extreme value copula of the PPPP copula
#'
#' Joint cdf of the upper extreme value copula of the PPPP copula
#' @param u Probability or first copula argument.
#' @param v Second copula argument.
#' @param be Finite positive parameter. If b <= be the copula is independence. The case a - b + be = 0 uses its analytic limit.
#' @param a Finite positive parameter. If b <= be the copula is independence. The case a - b + be = 0 uses its analytic limit.
#' @param b Finite positive parameter. If b <= be the copula is independence. The case a - b + be = 0 uses its analytic limit.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' pUEV_PPPP_COP(0.2, 0.4, 1.2, 1, 1)
#' pUEV_PPPP_COP(0.3, 0.4, 1, 1, 2)
pUEV_PPPP_COP <- function(u, v, be, a, b) {
  .check_UEV_PPPP(u, v, be, a, b)
  if (b <= be) return(u * v)
  if (u == 0 || v == 0) return(0)
  if (u == 1) return(v)
  if (v == 1) return(u)
  exp(-A_UEV_PPPP_COP(-log(u), -log(v), be, a, b))
}
#' Density function of the upper extreme value copula for the PPPP copula
#'
#' Density function of the upper extreme value copula for the PPPP copula
#' @param u1 First copula argument.
#' @param u2 Second copula argument.
#' @param be Finite positive parameter. If b <= be the copula is independence. The case a - b + be = 0 uses its analytic limit.
#' @param a Finite positive parameter. If b <= be the copula is independence. The case a - b + be = 0 uses its analytic limit.
#' @param b Finite positive parameter. If b <= be the copula is independence. The case a - b + be = 0 uses its analytic limit.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' dUEV_PPPP_COP(0.3, 0.4, 0.6, 1, 1)
#' dUEV_PPPP_COP(0.3, 0.4, 1, 1, 2)
dUEV_PPPP_COP <- function(u1, u2, be, a, b) {
  .check_UEV_PPPP(u1, u2, be, a, b)
  if (b <= be) return(1)
  if (u1 == 0 || u1 == 1 || u2 == 0 || u2 == 1)
    stop("Density arguments must be strictly in (0,1) when b > be.")
  w <- .UEV_PPPP_bundle(-log(u1), -log(u2), be, a, b)
  # c(u,v) = C(u,v)/(uv) * (L_1 L_2 - L_12).
  exp(-w$value - log(u1) - log(u2)) *
    (w$first_min * w$first_max + w$negative_mixed)
}
