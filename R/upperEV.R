# upperEV.R upper extreme value copula

#' Joint cdf of the upper extreme value copula
#'
#' Joint cdf of the upper extreme value copula of the full-range tail dependence copula
#' @param u, v: inputs for the extreme value copula
#' @param b: the parameter of the extreme value copula
#' @keywords upper extreme value copula
#' @export
#' @examples
#' pUEV_GGEE_COP(0.3, 0.4, 1.2)
pUEV_GGEE_COP <- function(u, v, b)
{
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

#' Density function of the upper extreme value copula
#'
#' Density of the upper extreme value copula of the full-range tail dependence copula
#' @param u, v: inputs for the extreme value copula
#' @param b: the parameter of the extreme value copula
#' @keywords upper extreme value copula
#' @export
#' @examples
#' dUEV_GGEE_COP(0.3, 0.4, 1.2)
dUEV_GGEE_COP <- function(u, v, b)
{
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




