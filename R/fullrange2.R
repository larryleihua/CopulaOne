# fullrange2.R full-range tail dependence on both upper and lower tails

###################### 
xmin <- 1e-10
xmax <- 1e+10
umin <- 1e-06
umax <- 1 - umin

#' Uniform score function
#'
#' Transform data to uniform scores, this function is from the R package CopulaModel.
#' @param data:   data to be transformed.
#' @keywords uniform score
#' @export
#' @examples
#' data  <- cbind(rnorm(10), rnorm(10))
#' uscore(data)
uscore <- function(data, aunif = -0.5)
{
  if (is.vector(data))
  {
    nr <- length(data)
    us <- ((1:nr) + aunif)/(nr + 1 + 2 * aunif)
    jj <- rank(data)
    out <- us[jj]
  } else
  {
    nc <- ncol(data)
    nr <- nrow(data)
    out <- matrix(0, nr, nc)
    us <- ((1:nr) + aunif)/(nr + 1 + 2 * aunif)
    for (j in 1:nc)
    {
      jj <- rank(data[, j])
      tem <- us[jj]
      out[, j] <- tem
    }
  }
  out
}

# gaussian quadrature copied from R package CopulaModel
gausslegendre <- function (nq)
{
  if (nq <= 0 || nq > 70) 
    stop("nq between 1 and 70")
  nq1 <- nq + 1
  out <- .C("gauleg", x1 = as.double(0), x2 = as.double(1), 
           as.integer(nq), xq = as.double(rep(0, nq1)), wq = as.double(rep(0, nq1)))
  list(nodes = out$xq[-1], weights = out$wq[-1])
}



#' CDF of univariate margins of the GGEE model
#'
#' CDF of univariate margins of the GGEE model
#' @param x:   data input.
#' @param a,b:   parameters.
#' @param method:   method for numerical integration: "GQ" for Gaussian Quadrature.
#' @keywords CDF
#' @export
#' @examples
#' pGGEE(2, 1, 2)
pGGEE <- function(x, a, b, method="GQ")
{
  if(method=="GQ")
  {
    
  }else
  {
    out <- tryCatch(1 - a/(a + b) * Re(hypergeo::hypergeo(1, b, a + b + 1, 1 - x, 
                                                          tol = 1e-06, maxiter = 10000)), error = function(err) FALSE, warning = function(err) FALSE)
    if (!is.logical(out) && is.finite(out))
    {
      return(out)
    } else
    {
      cat("Warning! NA returned!", "\n")
      return(NA)
    }
  }
}

#' Density of univariate margins of the GGEE model
#'
#' Density of univariate margins of the GGEE model
#' @param x:   data input.
#' @param a,b:   parameters.
#' @keywords Density
#' @export
#' @examples
#' dGGEE(3, 1, 2)
dGGEE <- function(x, a, b)
{
  out <- tryCatch((a * b)/(a + b)/(a + b + 1) * Re(hypergeo::hypergeo(2, b + 1, 
    a + b + 2, 1 - x, tol = 1e-06, maxiter = 10000)), error = function(err) FALSE, 
    warning = function(err) FALSE)
  if (!is.logical(out) && is.finite(out)) 
    return(out) else
    {
    cat("Warning! NA returned!", "\n")
    return(NA)
  }
}


intg_jdGGEE <- function(y, x1, x2, a, b)
{
  tem1 <- (x1*y+1-y)^(-2)
  tem2 <- (x2*y+1-y)^(-2)
  tem3 <- ((1-y)^(a+1))*(y^(b+1))
  tem1*tem2*tem3
}
intg_jdGGEE <- Vectorize(intg_jdGGEE, "y")
# plot(intg_jdGGEE(seq(0.0, 1, length=100), x1=2, x2=4, a=1.4, b=1.9))

#' Joint density of the GGEE model
#'
#' Joint density of the GGEE model
#' @param x1,x2:   data input.
#' @param a,b:   parameters.
#' @param flag: flag used in appell::appellf1()
#' @param integration   if T: use numerical integration instead of appellf1, default is F
#' @keywords Joint density
#' @export
#' @examples
#' jdGGEE(10,2, 1, 1)
#' jdGGEE(10,2, 1, 1, integration = T)
jdGGEE <- function(x1, x2, a, b, flag = 1, integration = F)
{
  if(integration==F)
  {
    tem1 <- a * b * (a + 1) * (b + 1)/(a + b)/(a + b + 1)/(a + b + 2)/(a + b + 3)
    tem2 <- tryCatch(Re(appell::appellf1(b + 2, 2, 2, a + b + 4, 1 - x1, 1 - x2, 
        userflag = flag)$val), error = function(err) FALSE, warning = function(err) FALSE)
    if (!is.logical(tem2) && is.finite(tem2)) 
    return(tem1 * tem2) else
    {
      cat("Warning! NA returned!", "\n")
      return(NA)
    }  
  }else{
    tmp <- tryCatch(integrate(intg_jdGGEE, lower = 0, upper = 1, x1 = x1, x2 = x2, 
      a = a, b = b, stop.on.error = T), error = function(err) FALSE, warning = function(err) FALSE)
    if (!is.logical(tmp))
    {
      intg <- tmp$value
      return(intg/beta(a, b))
    } else
    {
      cat("Warning! NA returned! (jdGGEE)", "\n")
      return(NA)
    }
  }
}

#' Quantile of univariate margins of the GGEE model
#'
#' Quantiel of univariate margins of the GGEE model
#' @param u:   quantile between 0 and 1.
#' @param a,b:   parameters.
#' @keywords Quantile
#' @export
#' @examples
#' qGGEE(0.9, 1, 1)
qGGEE <- function(u, a, b)
{
  if (u == 0)
  {
    out <- 0
  } else
  {
    tol <- 1e-06
    CDF <- -0.01
    DEN <- 1
    maxiter <- 1000
    kount <- 0
    t <- 0
    
    # -------------------------------- Now use modified Newton-Raphson
    # --------------------------------
    lower <- -1e+20
    upper <- 1e+20
    
    while ((kount < maxiter) && (abs(u - CDF) > tol))
    {
      kount <- kount + 1
      t <- t - (CDF - u)/DEN
      if (t < lower || t > upper)
      {
        t <- 0.5 * (lower + upper)
      }
      DEN <- dGGEE(exp(t), a, b) * exp(t)
      CDF <- pGGEE(exp(t), a, b)
      if (CDF < u)
      {
        lower <- t
      } else
      {
        upper <- t
      }
    }
    out <- exp(t)
  }
  return(out)
}

intg_Dx2 <- function(r, x1, x2, a, b)
{
  tem11 <- (x1 + r)^(-1)
  tem12 <- (x2 + r)^(-2)
  tem2 <- r^(a + 1)
  tem3 <- (1 + r)^(a + b)
  return(tem11 * tem12 * tem2/tem3)
}

Dx2 <- function(x1, x2, a, b)
{
  tem1 <- dGGEE(x2, a, b)
  if (!is.na(tem1))
  {
    intg_Dx2 <- Vectorize(intg_Dx2, "r")
    tmp <- tryCatch(integrate(intg_Dx2, lower = 0, upper = Inf, x1 = x1, x2 = x2, 
      a = a, b = b, stop.on.error = T), error = function(err) FALSE, warning = function(err) FALSE)
    if (!is.logical(tmp))
    {
      intg <- tmp$value
      return(tem1 - intg/beta(a, b))
    } else
    {
      cat("Warning! NA returned! (Dx2)", "\n")
      return(NA)
    }
  } else
  {
    cat("Warning! NA returned! (Dx2)", "\n")
    return(NA)
  }
}

#' Sampling based on the full-range tail dependence copula
#'
#' Generating random samples based on the bivariate copula that has full-range tail dependence in both upper and lower tails
#' @param n: sample size to be generated.
#' @param a,b: the two shape parameters.
#' @param seed: seed for randomness, and default is 1
#' @keywords simulation
#' @export
#' @examples
#' rGGEE_COP(20, 1.2, 0.2, seed = 100)
rGGEE_COP <- function(n, a, b, seed = 1)
{
  set.seed(seed)
  R1 <- rgamma(n, shape=a, 1)
  R2 <- rgamma(n, shape=b, 1)
  R <- R1 / R2
  H1 <- rexp(n)
  H2 <- rexp(n)
  H11 <- rexp(n)
  H22 <- rexp(n)
  X11 <- R * H1 / H11
  X22 <- R * H2 / H22
  
  u <- pGGEE(X11,a,b)
  v <- pGGEE(X22,a,b)
  cbind(u,v)
}

C2GGEE_COP <- function(u, v, a, b)
{
  x1 <- tryCatch(qGGEE(u, a, b), error = function(err) FALSE, warning = function(err) FALSE)
  x2 <- tryCatch(qGGEE(v, a, b), error = function(err) FALSE, warning = function(err) FALSE)
  if (is.logical(x1) || is.logical(x1))
  {
    return(NA)
    cat("Warning! NA returned! (C2GGEE_COP: qGGEE error!)", "\n")
  } else
  {
    if (is.finite(x1) && is.finite(x1))
    {
      tem1 <- Dx2(x1, x2, a, b)
      tem2 <- dGGEE(x2, a, b)
      if (is.finite(tem1) && is.finite(tem1))
      {
        return(tem1/tem2)
      } else
      {
        cat("Warning! NA returned! (C2GGEE_COP)", "\n")
        return(NA)
      }
    } else
    {
      cat("Warning! NA returned! (C2GGEE_COP)", "\n")
      return(NA)
    }
  }
}

dGGEE_COP_0 <- function(u, v, a, b, flag = 1, integration = F)
{
  q1 <- tryCatch(qGGEE(u, a, b), error = function(err) FALSE, warning = function(err) FALSE)
  q2 <- tryCatch(qGGEE(v, a, b), error = function(err) FALSE, warning = function(err) FALSE)
  if (is.logical(q1) || is.logical(q1))
  {
    return(NA)
    cat("Warning! NA returned! (dGGEE_COP: qGGEE error!)", "\n")
  } else
  {
    if (is.finite(q1) && is.finite(q2))
    {
      tem1 <- jdGGEE(q1, q2, a, b, flag, integration)
      tem2 <- dGGEE(q1, a, b)
      tem3 <- dGGEE(q2, a, b)
      if (is.finite(tem1) && is.finite(tem2) && is.finite(tem3))
      {
        return(tem1/tem2/tem3)
      } else
      {
        cat("Warning! NA returned! (dGGEE_COP)", "\n")
        return(NA)
      }
    } else
    {
      cat("Warning! NA returned! (dGGEE_COP)", "\n")
      return(NA)
    }
  }
}


#' Copula Density Function
#'
#' Copula density function of the bivariate copula that has full-range tail dependence in both upper and lower tails
#' @param u,v    values in (0,1).
#' @param a,b    the two shape parameters.
#' @param flag   used in the Appell's F1 function appellf1() of the R package 'appell'.
#' @param integration   if T: use numerical integration instead of appellf1, default is F
#' @keywords copula density
#' @export
#' @examples
#' dGGEE_COP(0.2, 0.4, 1.2, 0.2)
dGGEE_COP <- Vectorize(dGGEE_COP_0, c("u", "v"))

intg_jpGGEE <- function(y, x1, x2, a, b)
{
  tem1 <- (x1*y+1-y)^(-1)
  tem2 <- (x2*y+1-y)^(-1)
  tem3 <- ((1-y)^(a+1))*(y^(b-1))
  tem1*tem2*tem3
}
intg_jpGGEE <- Vectorize(intg_jpGGEE, "y")
# plot(intg_jpGGEE(seq(0.0, 1, length=100), x1=2, x2=4, a=1.4, b=1.9))

#' Joint CDF of the GGEE model
#'
#' Joint CDF of the GGEE model
#' @param x1,x2:   data input.
#' @param a,b:   parameters.
#' @keywords Joint CDF
#' @export
#' @examples
#' jpGGEE(0.2, 0.4, 1, 1)
#' jpGGEE(0.2, 0.4, 1, 1, integration = T)
jpGGEE <- function(x1,x2,a,b,flag=1, integration = F)
{
  tem1 <- pGGEE(x1,a,b)
  tem2 <- pGGEE(x2,a,b)
  if(integration==F)
  {
    tem3 <- a*(a+1)/(a+b)/(a+b+1)
    tem4 <- Re(appell::appellf1(b,1,1,a+b+2,1-x1,1-x2, userflag = flag)$val)
    out <- tem1 + tem2 - 1 + tem3 * tem4
    return(out)
  }else{
    tmp <- tryCatch(integrate(intg_jpGGEE, lower = 0, upper = 1, x1 = x1, x2 = x2, 
      a = a, b = b, stop.on.error = T), error = function(err) FALSE, warning = function(err) FALSE)
    if (!is.logical(tmp))
    {
      intg <- tmp$value
      return(tem1+tem2-1+intg/beta(a, b))
    }else
    {
      cat("Warning! NA returned! (jpGGEE)", "\n")
      return(NA)
    }    
  }
}



#' Joint CDF of the GGEE copula model
#'
#' Joint CDF of the GGEE copula model
#' @param u,v:   data input.
#' @param a,b:   parameters.
#' @param flag: flag used in the appell::appellf1() function.
#' @keywords Joint CDF
#' @export
#' @examples
#' pGGEE_COP(0.9, 0.3, 1, 1)
#' pGGEE_COP(0.9, 0.3, 1, 1, integration=T)
pGGEE_COP <- function(u,v,a,b,flag=1, integration=F)
{
  q1 <- qGGEE(u,a,b)
  q2 <- qGGEE(v,a,b)
  jpGGEE(q1,q2,a,b,flag = flag, integration=integration)
}

intg_tau_E <- function(xy,a,b)
{
  x <- xy[1]
  y <- xy[2]
  
  if(x==1 | y==1 | x==0 | y==0)
  {
   return(0) 
  }else
  {
    tem1 <- (-y*log((-1+x)*y/((-1+y)*x))+y*log((-1+x)*y/((-1+y)*x))*x-x+y)^2 
    tem2 <- x*(1-x)^(a-1)*(1-y)^(1+a)*(x*y)^b
    tem3 <- (y*(y^2-2*x*y+x^2)^2)
  }
  
  if(x == y)
  {
    out <- (1/4)*(((1-x)*(1-y))^(a-1))*((x*y)^(b-1))   
  }else
  {
    out <- tem1 * tem2 / tem3  
  }
  if(is.finite(out)) out else 0
}

#' Kendall's tau of the GGEE copula
#'
#' Kendall's tau of the bivariate copula that has full-range tail dependence in both upper and lower tails
#' @param a,b    the two shape parameters.
#' @keywords Kendall's tau
#' @export
#' @examples
#' tauGGEE_COP(0.5, 1)
tauGGEE_COP <- function(a,b,method=1)
{
  tmp <- R2Cuba::cuhre(2,1,intg_tau_E,a=a,b=b,lower = c(0,0), upper = c(1,1), flags=list(verbose=0))
  
  if(method==1)
  {
    tmp <- try(R2Cuba::cuhre(2,1,intg_tau_E,a=a,b=b,lower = c(0,0), upper = c(1,1), flags=list(verbose=0)), silent = T)
    if(is(tmp,"try-error"))
    {
      tmp2 <- try(cubature::adaptIntegrate(intg_tau_E, a=a, b=b, lowerLimit = c(0,0), upperLimit = c(1,1)), silent = T)
      if(is(tmp2,"try-error")){cat("tauGGEE_COP error at: ", a, b, "\n"); return(NA)}else{intg <- tmp2$integral}
    }else{intg <- tmp$value}
  }else if(method==2)
  {
    tmp <- cubature::adaptIntegrate(intg_tau_E, a=a, b=b, lowerLimit = c(0,0), upperLimit = c(1,1))
    intg <- tmp$integral
  }
  out <- 4 * intg/ (beta(a,b)^2) - 1
  out
}

intg_spr_C = function(u1u2,a,b,flag=1, integration=F)
{
  u1 <- u1u2[1]
  u2 <- u1u2[2]
  out <- 1 - u1 - u2 + pGGEE_COP(u1,u2,a,b,flag,integration)
  if(is.finite(out)) out else 0
}
#intg_spr_C <- Vectorize(intg_spr_C_0, "u1u2")

#' Spearman's rho of the GGEE copula
#'
#' Spearman's rho of the GGEE copula
#' @param a,b:    the two shape parameters.
#' @param flag: flag used in the appell::appellf1() function (default is -1L for auto choices)
#' @param method: different methods for the numerical integration. 1: R2Cuba::cuhre(); 2: cubature::adaptIntegrate(). 
#' @keywords Spearman's rho
#' @export
#' @examples
#' sprGGEE_COP(1.2, 0.6)
#' sprGGEE_COP(1.2, 0.6, method=2)
sprGGEE_COP <- function(a,b,flag=1, method=1, integration=F)
{
  if(method==1)
  {
    tmp <- try(R2Cuba::cuhre(2,1,intg_spr_C,a=a,b=b,flag=flag,integration=integration,lower = c(0,0), upper = c(1,1), flags=list(verbose=0)), silent = T)
    if(is(tmp,"try-error"))
    {
      tmp2 <- try(cubature::adaptIntegrate(intg_spr_C, a=a, b=b, flag=flag, integration=integration,lowerLimit = c(0,0), upperLimit = c(1,1)), silent = T)
      if(is(tmp2,"try-error")){cat("sprGGEE_COP error at: ", a, b, "\n"); return(NA)}else{intg <- tmp2$integral}
    }else{intg <- tmp$value}
  }else if(method==2)
  {
    tmp <- cubature::adaptIntegrate(intg_spr_C, a=a, b=b, flag=flag,integration=integration,lowerLimit = c(0,0), upperLimit = c(1,1))
    intg <- tmp$integral
  }
  out <- 12 * intg - 3
  out
}
