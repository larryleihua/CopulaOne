# fullrange2.R full-range tail dependence on both upper and lower tails

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

#' CDF of univariate margins of the GGEE model
#'
#' CDF of univariate margins of the GGEE model
#' @param x:   data input.
#' @param al,be:   parameters.
#' @keywords CDF
#' @export
#' @examples
#' pGGEE(2, 1, 2)
pGGEE <- function(x, al, be, maxit=100000)
{
  idx1 <- (x >= 1.0e2)
  idx2 <- (x <= 1.0e2)
  hyp1 <- Re(hypergeo::hypergeo(1, be, al + be + 1, 1 - x[idx1], tol = 1e-06, maxiter = maxit))
  hyp2 <- Re(hypergeo::hypergeo(1, al+1, al + be + 1, 1 - 1/x[idx2], tol = 1e-06, maxiter = maxit))/x[idx2] # Euler and Pfaff transformations
  hyp <- rep(0, length(x))
  hyp[idx1] <- hyp1
  hyp[idx2] <- hyp2
  out <- tryCatch(1 - al/(al + be) * hyp, error = function(err) FALSE, warning = function(err) FALSE)
  if (all(!is.logical(out) & is.finite(out)))
  {
    return(out)
  } else
  {
    cat("Warning! NA returned!", "\n")
    return(NA)
  }
}

#' Density of univariate margins of the GGEE model
#' 
#' Density of univariate margins of the GGEE model
#' @param x:   data input.
#' @param al,be:   parameters.
#' @keywords Density
#' @export
#' @examples
#' dGGEE(3, 1, 2)
dGGEE <- function(x, al, be, maxit=100000)
{
  out <- tryCatch((al * be)/(al + be)/(al + be + 1) * Re(hypergeo::hypergeo(2, be + 1, 
                                                                            al + be + 2, 1 - x, tol = 1e-06, maxiter = maxit)), error = function(err) FALSE, 
                  warning = function(err) FALSE)
  if (all(!is.logical(out) & is.finite(out)))
    return(out) else
    {
      cat("Warning! NA returned!", "\n")
      return(NA)
    }
}

intg_jdGGEE <- function(y, x1, x2, al, be)
{
  tem1 <- (x1*y+1-y)^(-2)
  tem2 <- (x2*y+1-y)^(-2)
  tem3 <- ((1-y)^(al+1))*(y^(be+1))
  tem1*tem2*tem3
}
intg_jdGGEE <- Vectorize(intg_jdGGEE, "y")


#' Joint density of the GGEE model
#' 
#' Joint density of the GGEE model
#' @param x1,x2:   data input.
#' @param al,be:   parameters.
#' @param flag: flag used in appell::appellf1()
#' @param integration   if T: use numerical integration instead of appellf1, default is F
#' @keywords Joint density
#' @export
#' @examples
#' jdGGEE(10,2, 1, 1)
#' jdGGEE(10,2, 1, 1, integration = T)
jdGGEE <- function(x1, x2, al, be, flag = 1, integration = F)
{
  if(integration==F)
  {
    tem1 <- al * be * (al + 1) * (be + 1)/(al + be)/(al + be + 1)/(al + be + 2)/(al + be + 3)
    tem2 <- tryCatch(Re(appell::appellf1(be + 2, 2, 2, al + be + 4, 1 - x1, 1 - x2, 
                                         userflag = flag)$val), error = function(err) FALSE, warning = function(err) FALSE)
    if (all(!is.logical(tem2) & is.finite(tem2)))
      return(tem1 * tem2) else
      {
        cat("Warning! NA returned!", "\n")
        return(NA)
      }  
  }else{
    tmp <- tryCatch(integrate(intg_jdGGEE, lower = 0, upper = 1, x1 = x1, x2 = x2, 
                              al = al, be = be, stop.on.error = T), error = function(err) FALSE, warning = function(err) FALSE)
    if (!is.logical(tmp))
    {
      intg <- tmp$value
      return(intg/beta(al, be))
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
#' @param al,be:   parameters.
#' @keywords Quantile
#' @export
#' @examples
#' qGGEE(0.9, 1, 1)
qGGEE <- function(u, al, be, maxit=100000)
{
  if (u == 0)
  {
    out <- 0
  } else
  {
    tol <- 1e-06
    CDF <- -0.01
    DEN <- 1
    kount <- 0
    t <- 0
    
    # -------------------------------- Now use modified Newton-Raphson
    # --------------------------------
    lower <- -1e+20
    upper <- 1e+20
    
    while ((kount < maxit) && (abs(u - CDF) > tol))
    {
      kount <- kount + 1
      t <- t - (CDF - u)/DEN
      if (t < lower || t > upper)
      {
        t <- 0.5 * (lower + upper)
      }
      DEN <- dGGEE(exp(t), al, be, maxit=maxit) * exp(t)
      CDF <- pGGEE(exp(t), al, be, maxit=maxit)
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

intg_Dx2_GGEE <- function(r, x1, x2, al, be)
{
  tem11 <- (x1 + r)^(-1)
  tem12 <- (x2 + r)^(-2)
  tem2 <- r^(al + 1)
  tem3 <- (1 + r)^(al + be)
  return(tem11 * tem12 * tem2/tem3)
}

Dx2_GGEE <- function(x1, x2, al, be)
{
  tem1 <- dGGEE(x2, al, be)
  if (!is.na(tem1))
  {
    intg_Dx2_GGEE <- Vectorize(intg_Dx2_GGEE, "r")
    tmp <- tryCatch(integrate(intg_Dx2_GGEE, lower = 0, upper = Inf, x1 = x1, x2 = x2, 
                              al = al, be = be, stop.on.error = T), error = function(err) FALSE, warning = function(err) FALSE)
    if (!is.logical(tmp))
    {
      intg <- tmp$value
      return(tem1 - intg/beta(al, be))
    } else
    {
      cat("Warning! NA returned! (Dx2_GGEE)", "\n")
      return(NA)
    }
  } else
  {
    cat("Warning! NA returned! (Dx2_GGEE)", "\n")
    return(NA)
  }
}

#' Sampling based on the full-range tail dependence copula
#'
#' Generating random samples based on the bivariate copula that has full-range tail dependence in both upper and lower tails
#' @param n: sample size to be generated.
#' @param al,be: the two shape parameters.
#' @param seed: seed for randomness
#' @keywords simulation
#' @export
#' @examples
#' rGGEE_COP(20, 1.2, 0.2, seed = 100)
rGGEE_COP <- function(n, al, be, seed = NULL, maxit=100000)
{
  if(!is.null(seed)){ set.seed(seed) }
  R1 <- rgamma(n, shape=al, 1) # rgamma can generate 0
  R2 <- rgamma(n, shape=be, 1)
  R1[R1==0] <- .Machine$double.xmin
  R2[R2==0] <- .Machine$double.xmin
  R <- R1 / R2
  H1 <- rexp(n)
  H2 <- rexp(n)
  H11 <- rexp(n)
  H22 <- rexp(n)
  X11 <- R * H1 / H11
  X22 <- R * H2 / H22
  X11[X11==0] <- .Machine$double.xmin
  X22[X22==0] <- .Machine$double.xmin
  X11[X11==Inf] <- .Machine$double.xmax
  X22[X22==Inf] <- .Machine$double.xmax
  
  u <- pGGEE(X11,al,be, maxit=maxit)
  v <- pGGEE(X22,al,be, maxit=maxit)
  cbind(u,v)
}

#' Conditional GGEE copula
#'
#' Partial derivative wrt the second argument of the GGEE copula
#' @param u,v:   data input.
#' @param al,be:   parameters.
#' @keywords conditional copula
#' @export
#' @examples
#' C2GGEE_COP(0.2, 0.6, al = 1.2, be = 0.8)
C2GGEE_COP <- function(u, v, al, be, maxit=100000)
{
  x1 <- tryCatch(qGGEE(u, al, be, maxit=maxit), error = function(err) FALSE, warning = function(err) FALSE)
  x2 <- tryCatch(qGGEE(v, al, be, maxit=maxit), error = function(err) FALSE, warning = function(err) FALSE)
  if (is.logical(x1) || is.logical(x2))
  {
    return(NA)
    cat("Warning! NA returned! (C2GGEE_COP: qGGEE error!)", "\n")
  } else
  {
    if (all(is.finite(x1) & is.finite(x2)))
    {
      tem1 <- Dx2_GGEE(x1, x2, al, be)
      tem2 <- dGGEE(x2, al, be, maxit=maxit)
      if (all(is.finite(tem1) & is.finite(tem2)))
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

dGGEE_COP_0 <- function(u, v, al, be, flag = 1, integration = F, maxit=100000)
{
  q1 <- tryCatch(qGGEE(u, al, be, maxit=maxit), error = function(err) FALSE, warning = function(err) FALSE)
  q2 <- tryCatch(qGGEE(v, al, be, maxit=maxit), error = function(err) FALSE, warning = function(err) FALSE)
  if (is.logical(q1) || is.logical(q1))
  {
    return(NA)
    cat("Warning! NA returned! (dGGEE_COP: qGGEE error!)", "\n")
  } else
  {
    if (all(is.finite(q1) & is.finite(q2)))
    {
      tem1 <- jdGGEE(q1, q2, al, be, flag, integration)
      tem2 <- dGGEE(q1, al, be, maxit=maxit)
      tem3 <- dGGEE(q2, al, be, maxit=maxit)
      if (all(is.finite(tem1) & is.finite(tem2) & is.finite(tem3)))
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


#' Copula Density Function - GGEE_COP
#'
#' Copula density function of the bivariate copula that has full-range tail dependence in both upper and lower tails
#' @param u,v    values in (0,1).
#' @param al,be    the two shape parameters.
#' @param flag   used in the Appell's F1 function appellf1() of the R package 'appell'.
#' @param integration   if T: use numerical integration instead of appellf1, default is F
#' @keywords copula density
#' @export
#' @examples
#' dGGEE_COP(0.2, 0.4, 1.2, 0.2)
dGGEE_COP <- Vectorize(dGGEE_COP_0, c("u", "v"))

intg_jpGGEE <- function(y, x1, x2, al, be)
{
  tem1 <- (x1*y+1-y)^(-1)
  tem2 <- (x2*y+1-y)^(-1)
  tem3 <- ((1-y)^(al+1))*(y^(be-1))
  tem1*tem2*tem3
}
intg_jpGGEE <- Vectorize(intg_jpGGEE, "y")
# plot(intg_jpGGEE(seq(0.0, 1, length=100), x1=2, x2=4, al=1.4, be=1.9))

#' Joint CDF of the GGEE model
#'
#' Joint CDF of the GGEE model
#' @param x1,x2:   data input.
#' @param al,be:   parameters.
#' @keywords Joint CDF
#' @export
#' @examples
#' jpGGEE(0.2, 0.4, 1, 1)
#' jpGGEE(0.2, 0.4, 1, 1, integration = T)
jpGGEE <- function(x1,x2,al,be,flag=1, integration = F, maxit=100000)
{
  tem1 <- pGGEE(x1,al,be, maxit=maxit)
  tem2 <- pGGEE(x2,al,be, maxit=maxit)
  if(integration==F)
  {
    tem3 <- al*(al+1)/(al+be)/(al+be+1)
    tem4 <- Re(appell::appellf1(be,1,1,al+be+2,1-x1,1-x2, userflag = flag)$val)
    out <- tem1 + tem2 - 1 + tem3 * tem4
    return(out)
  }else{
    tmp <- tryCatch(integrate(intg_jpGGEE, lower = 0, upper = 1, x1 = x1, x2 = x2, 
                              al = al, be = be, stop.on.error = T), error = function(err) FALSE, warning = function(err) FALSE)
    if (!is.logical(tmp))
    {
      intg <- tmp$value
      return(tem1+tem2-1+intg/beta(al, be))
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
#' @param al,be:   parameters.
#' @param flag: flag used in the appell::appellf1() function.
#' @keywords Joint CDF
#' @export
#' @examples
#' pGGEE_COP(0.9, 0.3, 1, 1)
#' pGGEE_COP(0.9, 0.3, 1, 1, integration=T)
pGGEE_COP <- function(u,v,al,be,flag=1, integration=F, maxit=100000)
{
  q1 <- qGGEE(u,al,be,maxit=maxit)
  q2 <- qGGEE(v,al,be,maxit=maxit)
  jpGGEE(q1,q2,al,be,flag = flag, integration=integration, maxit=maxit)
}

intg_tau_E <- function(xy,al,be)
{
  x <- xy[1]
  y <- xy[2]
  
  if(x==1 | y==1 | x==0 | y==0)
  {
    return(0) 
  }else
  {
    tem1 <- (-y*log((-1+x)*y/((-1+y)*x))+y*log((-1+x)*y/((-1+y)*x))*x-x+y)^2 
    tem2 <- x*(1-x)^(al-1)*(1-y)^(1+al)*(x*y)^be
    tem3 <- (y*(y^2-2*x*y+x^2)^2)
  }
  
  if(x == y)
  {
    out <- (1/4)*(((1-x)*(1-y))^(al-1))*((x*y)^(be-1))   
  }else
  {
    out <- tem1 * tem2 / tem3  
  }
  if(is.finite(out)) out else 0
}

#' Kendall's tau of the GGEE copula
#'
#' Kendall's tau of the bivariate copula that has full-range tail dependence in both upper and lower tails
#' @param al,be    the two shape parameters.
#' @keywords Kendall's tau
#' @export
#' @examples
#' tauGGEE_COP(0.5, 1)
tauGGEE_COP <- function(al,be)
{
   cat("WARNING: please use tauGGEE_COP_sim for smaller values (< 0.2) of al and be!!", "\n")
   tmp <- try(cubature::adaptIntegrate(intg_tau_E, al=al, be=be, lowerLimit = c(0,0), upperLimit = c(1,1)), silent = T)
   if(is(tmp,"try-error"))
   {
     cat("tauGGEE_COP error at: ", al, be, "\n"); return(NA)
   }else{intg <- tmp$integral}
  out <- 4 * intg/ (beta(al,be)^2) - 1
  out
}

#' Kendall's tau of the GGEE copula (based on simulations)
#'
#' Kendall's tau of the bivariate copula that has full-range tail dependence in both upper and lower tails
#' @param al,be    the two shape parameters.
#' @keywords Kendall's tau
#' @export
#' @examples
#' tauGGEE_COP_sim(0.5, 1)
tauGGEE_COP_sim <- function(al, be, n=10000000)
{
  R1 <- rgamma(n, shape=al, 1)
  R2 <- rgamma(n, shape=be, 1)
  R1[R1==0] <- .Machine$double.xmin
  R2[R2==0] <- .Machine$double.xmin
  R <- R1 / R2
  H1 <- rexp(n)
  H2 <- rexp(n)
  H11 <- rexp(n)
  H22 <- rexp(n)
  X11 <- R * H1 / H11
  X22 <- R * H2 / H22
  X11[X11==0] <- .Machine$double.xmin
  X22[X22==0] <- .Machine$double.xmin
  X11[X11==Inf] <- .Machine$double.xmax
  X22[X22==Inf] <- .Machine$double.xmax
  
  R1_ <- rgamma(n, shape=al, 1)
  R2_ <- rgamma(n, shape=be, 1)
  R1_[R1_==0] <- .Machine$double.xmin
  R2_[R2_==0] <- .Machine$double.xmin
  R_ <- R1_ / R2_
  H1_ <- rexp(n)
  H2_ <- rexp(n)
  H11_ <- rexp(n)
  H22_ <- rexp(n)
  X11_ <- R_ * H1_ / H11_
  X22_ <- R_ * H2_ / H22_
  X11_[X11_==0] <- .Machine$double.xmin
  X22_[X22_==0] <- .Machine$double.xmin
  X11_[X11_==Inf] <- .Machine$double.xmax
  X22_[X22_==Inf] <- .Machine$double.xmax
  
  out <- 4*(sum((X11 <= X11_)&(X22 <= X22_)) / n) - 1
  out
}


intg_spr_C = function(u1u2,al,be,flag=1, integration=F, maxit=100000)
{
  u1 <- u1u2[1]
  u2 <- u1u2[2]
  out <- 1 - u1 - u2 + pGGEE_COP(u1,u2,al,be,flag,integration, maxit=maxit)
  if(is.finite(out)) out else 0
}
#intg_spr_C <- Vectorize(intg_spr_C_0, "u1u2")

#' Spearman's rho of the GGEE copula
#'
#' Spearman's rho of the GGEE copula
#' @param al,be:    the two shape parameters.
#' @param flag: flag used in the appell::appellf1() function (default is -1L for auto choices)
#' @keywords Spearman's rho
#' @export
#' @examples
#' sprGGEE_COP(1.2, 0.6)
sprGGEE_COP <- function(al,be,flag=1, integration=F)
{
  tmp <- try(cubature::adaptIntegrate(intg_spr_C, al=al, be=be, flag=flag,integration=integration,lowerLimit = c(0,0), upperLimit = c(1,1)), silent = T)
  if(is(tmp,"try-error"))
  {
   cat("sprGGEE_COP error at: ", al, be, "\n"); return(NA)    
  }else{intg <- tmp$integral}
  out <- 12 * intg - 3
  out
}

############################################################
## full-range tail dependence copula with Pareto mixtures ##
############################################################

HX_PPPP <- function(mu, th1, th2, ga1, ga2)
{
  tem1 <- ga1*ga2
  tem2 <- ga1+ga2
  tem3 <- mu+ga1
  tem4 <- mu-ga2
  
  if(mu!=-ga1 & mu!=ga2)
  {
    if( 1 >= th2 )
    {
      out <- tem1/tem2/tem3*(th2^tem3 - th1^tem3)
    }else if( 1 < th2 & 1>= th1)
    {
      out <- tem1/tem2*((th2^tem4)/tem4 - (th1^tem3)/tem3) - tem1/tem3/tem4
    }else # 1 < th1
    {
      out <- tem1/tem2/tem4*(th2^tem4 - th1^tem4)
    }
  }else if(mu == -ga1)
  {
    if( 1 >= th2 )
    {
      out <- tem1/tem2*(log(th2) - log(th1))
    }else if( 1 < th2 & 1>= th1)
    {
      out <- tem1/tem2*((th2^tem4-1)/tem4 - log(th1))
    }else # 1 < th1
    {
      out <- tem1/tem2/tem4*(th2^tem4 - th1^tem4)
    }
  }else # mu == ga2
  {
    if( 1 >= th2 )
    {
      out <- tem1/tem2/tem3*(th2^tem3 - th1^tem3)
    }else if( 1 < th2 & 1>= th1)
    {
      out <- tem1/tem2*(log(th2) - (th1^tem3-1)/tem3)
    }else # 1 < th1
    {
      out <- tem1/tem2*(log(th2) - log(th1))
    }
  }
  return(out)
}
# HX_PPPP(2, 0.6, 1.7, 1.2, 0.7)

#' CDF of univariate margins of the PP model (i.e., Pareto I / Pareto I)
#' 
#' CDF of univariate margins of the PP model (i.e., Pareto I / Pareto I)
#' @param x:   data input.
#' @param al,be:   parameters.
#' @keywords CDF
#' @export
#' @examples
#' plot(sapply(seq(0,10, length=50), function(x){pPP(x, 2,2)}), type="l")
pPP <- function(x, ga1, ga2)
{
  if(x>=0 & x<1)
  {
    out <- ga2/(ga1+ga2)*(x^ga1)
  }else
  {
    out <- 1 - ga1/(ga1+ga2)*(x^(-ga2))
  }
  return(out)
}

#' density of univariate margins of the PP model (i.e., Pareto I / Pareto I)
#' 
#' density of univariate margins of the PP model (i.e., Pareto I / Pareto I)
#' @param x:   data input.
#' @param al,be:   parameters.
#' @keywords pdf
#' @export
#' @examples
#' plot(sapply(seq(0,10, length=50), function(x){dPP(x, 2,2)}), type="l")
dPP <- function(x, ga1, ga2)
{
  if(x>=0 & x<1)
  {
    out <- ga1*ga2/(ga1+ga2)*(x^(ga1-1))
  }else
  {
    out <- ga1*ga2/(ga1+ga2)*(x^(-ga2-1))
  }
  return(out)
}

#' quantile of univariate margins of the PP model (i.e., Pareto I / Pareto I)
#' 
#' quantile of univariate margins of the PP model (i.e., Pareto I / Pareto I)
#' @param x:   data input.
#' @param al,be:   parameters.
#' @keywords quantile
#' @export
#' @examples
#' plot(sapply(seq(0.99,1-0.001, length=1000), function(x){qPP(x, 0.2,0.3)}), type="l")
qPP <- function(u, ga1, ga2)
{
  out <- uniroot(function(x){pPP(x,ga1,ga2)-u}, c(0,9e99))$root
  return(out)
}

# internal function the g function for pPPPP()
g_PPPP <- function(a1,a2,a3,a4)
{
  1/((1+a1/a2)*(1+a1/a4)*(1-a1/a3))
}

#' CDF of univariate margins of the PPPP model
#' 
#' CDF of univariate margins of the PPPP model
#' @param x:   data input.
#' @param al,be,a,b:   parameters.
#' @param log: return log or not
#' @keywords CDF
#' @export
#' @examples
#' plot(sapply(seq(0.0001,20, length=50), function(x){pPPPP(x, 1,1,2,2,log=F)}), type="l",ylab="")
pPPPP <- function(x, al, be, a, b, log=F)
{
  if(x<=.Machine$double.xmin)
  {
    return(.Machine$double.xmin)
  }else
  {
    FR <- pPP(x,al,be)
    I1 <- HX_PPPP(b,0,x,al,be)
    I2 <- HX_PPPP(-a,x,Inf,al,be)
    out <- FR-a/(a+b)*(x^(-b))*I1+b/(a+b)*(x^(a))*I2
    if(log==F)
    {
      return(out)  
    }else
    {
      return(log(out))
    }
  }
}

# plot(sapply(seq(0,20, length=50), function(x){pPPPP_this_is_the_case_a_not_al_b_not_be(x, 1,1,2,2)}), type="l",ylab="")
pPPPP_this_is_the_case_a_not_al_b_not_be <- function(x, al, be, a, b)
{
  if(x>=0 & x<1)
  {
    out <- g_PPPP(al,be,a,b)*(x^al) + g_PPPP(a,b,al,be)*(x^a)
  }else
  {
    out <- 1 - g_PPPP(be,al,b,a)*(x^(-be)) - g_PPPP(b,a,be,al)*(x^(-b))
  }
  return(out)
}

#' quantile of univariate margins of the PPPP model
#' 
#' quantile of univariate margins of the PPPP model
#' @param u:   data input.
#' @param al,be,a,b:   parameters.
#' @param log: if use log in it but the result is still the original scale
#' @keywords quantile
#' @export
#' @examples
#' plot(sapply(seq(0.001,0.999, length=100), function(x){qPPPP(x, 0.3, 1.3, 1, 1,log=F)}), type="l",ylab="")

# try Newton
qPPPP <- function(u, al, be, a, b, maxit=100000)
{
  if (u == 0)
  {
    out <- 0
  } else
  {
    tol <- 1e-06
    CDF <- -0.01
    DEN <- 1
    maxiter <- maxit
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
      DEN <- dPPPP(exp(t), al, be, a, b) * exp(t)
      CDF <- pPPPP(exp(t), al, be, a, b)
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

# this is not reliable 
# qPPPP_not_reliable <- function(u, al, be, a, b, log=T)
# {
#   if(log==F)
#   {
#     out <- uniroot(function(x){pPPPP(x,al,be,a,b)-u}, c(0,9e99))$root  
#   }else
#   {
#     out <- uniroot(function(x){pPPPP(x,al,be,a,b,log=log)-log(u)}, c(1e-20,9e99))$root
#   }
#   return(out)
# }

#' Density of univariate margins of the PPPP model (a!=al, b!=be)
#' 
#' Density of univariate margins of the PPPP model (a!=al, b!=be)
#' @param x:   data input.
#' @param al,be,a,b:   parameters.
#' @keywords Density
#' @export
#' @examples
#' plot(sapply(seq(0,20, length=100), function(x){dPPPP(x, 1,1,2,2)}), type="l",ylab="")
dPPPP <- function(x, al, be, a, b)
{
  w <- a*b/(a+b)
  HR1 <- HX_PPPP(b,0,x,al,be)
  HR2 <- HX_PPPP(-a,x,Inf,al,be)
  out <- w*((x^(-b-1))*HR1+(x^(a-1))*HR2)
  return(out)
}

# plot(sapply(seq(0,20, length=100), function(x){dPPPP_this_is_the_case_a_not_al_b_not_be(x, 1,1,2,2)}), type="l",ylab="")
dPPPP_this_is_the_case_a_not_al_b_not_be <- function(x, al,be,a,b)
{
  if(x>=0 & x<1)
  {
    out <- al*g_PPPP(al,be,a,b)*(x^(al-1)) + a*g_PPPP(a,b,al,be)*(x^(a-1))
  }else
  {
    out <- be*g_PPPP(be,al,b,a)*(x^(-be-1)) + b*g_PPPP(b,a,be,al)*(x^(-b-1))
  }
  return(out)
}

#' Joint density of the PPPP model
#' 
#' Joint density of the PPPP model
#' @param x1,x2:   data input.
#' @param al,be,a,b:   parameters.
#' @keywords Joint density
#' @export
#' @examples
#' jdPPPP(10,2, 1, 1,2,2)
jdPPPP <- function(x1, x2, al,be,a,b)
{
  w <- a*b/(a+b)
  xm <- min(x1,x2)
  xp <- max(x1,x2)
  tem1 <- (xm*xp)^(-b-1)
  tem2 <- (xm^(a-1))*(xp^(-b-1))
  tem3 <- (xm*xp)^(a-1)
  HR1 <- HX_PPPP(2*b,0,xm,al,be)
  HR2 <- HX_PPPP(b-a,xm,xp,al,be)
  HR3 <- HX_PPPP(-2*a,xp,Inf,al,be)
  out <- (w^2)*(tem1*HR1+tem2*HR2+tem3*HR3)
  return(out)
}

dPPPP_COP_0 <- function(u,v,al,be,a,b)
{
  q1 <- tryCatch(qPPPP(u,al,be,a,b), error = function(err) FALSE, warning = function(err) FALSE)
  q2 <- tryCatch(qPPPP(v,al,be,a,b), error = function(err) FALSE, warning = function(err) FALSE)
  if (is.logical(q1) || is.logical(q1))
  {
    return(NA)
    cat("Warning! NA returned! (dPPPP_COP: qPPPP error!)", "\n")
  } else
  {
    if (all(is.finite(q1) & is.finite(q2)))
    {
      tem1 <- jdPPPP(q1,q2,al,be,a,b)
      tem2 <- dPPPP(q1,al,be,a,b)
      tem3 <- dPPPP(q2,al,be,a,b)
      if (all(is.finite(tem1) & is.finite(tem2) & is.finite(tem3)))
      {
        return(tem1/tem2/tem3)
      } else
      {
        cat("Warning! NA returned! (dPPPP_COP)", "\n")
        return(NA)
      }
    } else
    {
      cat("Warning! NA returned! (dPPPP_COP)", "\n")
      return(NA)
    }
  }
}

#' Copula Density Function - PPPP_COP
#'
#' PPPP_COP - Copula density function of the bivariate copula that has full-range tail dependence in both upper and lower tails
#' @param u,v    values in (0,1).
#' @param al,be,a,b    the four shape parameters.
#' @keywords copula density
#' @export
#' @examples
#' dPPPP_COP(0.2, 0.4, 1,1,2,2)
dPPPP_COP <- Vectorize(dPPPP_COP_0, c("u", "v"))

#' Copula Density Function of PPPP_COP (when a=b=1)
#'
#' PPPP_COP - Copula density function of PPPP_COP (when a=b=1)
#' @param u,v    values in (0,1).
#' @param al,be    the four shape parameters.
#' @export
#' @examples
#' dPPPP_COP_1(0.2, 0.4, 1,1)
#' dPPPP_COP_1(0.8, 0.3, 1,1)
#' dPPPP_COP_1(c(0.2,0.8), c(0.4, 0.3), 1,1)
dPPPP_COP_1 <- function(uvec, vvec, al, be)
{
  dPPPP_COP(uvec, vvec, al, be,1,1)
}

#' Copula Density Function of PPPP_COP that is rotated for 90 degrees counter clock-wise (when a=b=1)
#'
#' Copula density function of PPPP_COP that is rotated for 90 degrees counter clock-wise (when a=b=1)
#' @param u,v    values in (0,1).
#' @param al,be    the four shape parameters.
#' @export
#' @examples
#' dPPPP_COP_1_90(0.2, 0.4, 1,1)
#' dPPPP_COP_1_90(0.8, 0.3, 1,1)
#' dPPPP_COP_1_90(c(0.2,0.8), c(0.4, 0.3), 1,1)
dPPPP_COP_1_90 <- function(uvec, vvec, al, be)
{
  dPPPP_COP(vvec, 1-uvec, al, be,1,1)
}

# X1onX2_PPPP: P[X1<=x1|X2=x2], note that this is different than Dx2_PPPP, Dx2_PPPP = X1onX2_PPPP*dPPPP(X2)
X1onX2_PPPP <- function(x1, x2, al, be, a, b)
{
  if(x1 == 0){return(0)}else
  {
    if(x2 <= x1)
    {
      tem1 <- a*b*(x2^(-b-1))/(a+b)
      tem2 <- (a^2)*b*(x1^(-b))*(x2^(-b-1))/((a+b)^2)
      tem3 <- a*b*(x2^(a-1))/(a+b)
      tem4 <- (a^2)*b*(x1^(-b))*(x2^(a-1))/((a+b)^2)
      tem5 <- a*(b^2)*(x1^a)*(x2^(a-1))/((a+b)^2)
      HR1 <- HX_PPPP(b, 0, x2, al, be)
      HR2 <- HX_PPPP(2*b, 0, x2, al, be)
      HR3 <- HX_PPPP(-a, x2, x1, al, be)
      HR4 <- HX_PPPP(b-a, x2, x1, al, be)
      HR5 <- HX_PPPP(-2*a, x1, Inf, al, be)
      out <- (tem1*HR1 - tem2*HR2 + tem3*HR3 - tem4*HR4 + tem5*HR5) / dPPPP(x2, al, be, a, b)
    }else
    {
      tem1 <- a*b*(x2^(-b-1))/(a+b)
      tem2 <- (a^2)*b*(x1^(-b))*(x2^(-b-1))/((a+b)^2)
      tem3 <- a*(b^2)*(x1^a)*(x2^(-b-1))/((a+b)^2)
      tem4 <- a*(b^2)*(x1^a)*(x2^(a-1))/((a+b)^2)
      HR1 <- HX_PPPP(b, 0, x1, al, be)
      HR2 <- HX_PPPP(2*b, 0, x1, al, be)
      HR3 <- HX_PPPP(b-a, x1, x2, al, be)
      HR4 <- HX_PPPP(-2*a, x2, Inf, al, be)
      out <- (tem1*HR1 - tem2*HR2 + tem3*HR3 + tem4*HR4) / dPPPP(x2, al, be, a, b)
    }
    return(out)
  }
}
# plot(sapply(seq(0,20, length=50), function(x1){X1onX2_PPPP(x1, 2, 1,1,2,2)}), type="l")

#' Conditional PPPP copula
#'
#' Partial derivative wrt the second argument of the PPPP copula
#' @param u,v:   data input.
#' @param al,be,a,b:   parameters.
#' @keywords conditional copula
#' @export
#' @examples
#' C2PPPP_COP(0.2, 0.6, al=1.2, be=0.8, a=1, b=1)
C2PPPP_COP <- function(u, v, al, be, a, b)
{
  x1 <- tryCatch(qPPPP(u, al, be, a, b), error = function(err) FALSE, warning = function(err) FALSE)
  x2 <- tryCatch(qPPPP(v, al, be, a, b), error = function(err) FALSE, warning = function(err) FALSE)
  if (is.logical(x1) || is.logical(x2))
  {
    return(NA)
    cat("Warning! NA returned! (C2PPPP_COP: qPPPP error!)", "\n")
  } else
  {
    if (all(is.finite(x1) & is.finite(x2)))
    {
      tem1 <- X1onX2_PPPP(x1, x2, al, be, a, b)

      if (is.finite(tem1))
      {
        return(tem1)
      } else
      {
        cat("Warning! NA returned! (C2PPPP_COP)", "\n")
        return(NA)
      }
    } else
    {
      cat("Warning! NA returned! (C2PPPP_COP)", "\n")
      return(NA)
    }
  }
}
# plot(sapply(seq(0,1, length=50), function(u){C2PPPP_COP(u, 0.1, 0.5,1.5,1,1)}), type="l")

# quantile function of Pareto-I
qParetoI <- function(u, a)
{
  (1-u)^(-1/a)
}

#' Sampling based on the PPPP copula
#'
#' Generating random samples based on the PPPP copula
#' @param n: sample size to be generated.
#' @param al, be, a,b: the two shape parameters.
#' @param seed: seed for randomness
#' @keywords simulation
#' @export
#' @examples
#' plot(rPPPP_COP(100, 1.2, 0.2, 1, 1, seed = 100))
rPPPP_COP <- function(n, al, be, a, b, seed = NULL)
{
  if(!is.null(seed)){ set.seed(seed) }
  u1 <- runif(n)
  u2 <- runif(n)
  u3 <- runif(n)
  u4 <- runif(n)
  u5 <- runif(n)
  u6 <- runif(n)

  R1 <- sapply(u1, function(u){qParetoI(u, be)})
  R2 <- sapply(u2, function(u){qParetoI(u, al)})
  R1[R1==Inf] <- .Machine$double.xmax
  R2[R2==Inf] <- .Machine$double.xmax

  R <- R1 / R2
  
  H1 <- sapply(u3, function(u){qParetoI(u, b)})
  H11 <- sapply(u4, function(u){qParetoI(u, a)})
  H1[H1==Inf] <- .Machine$double.xmax
  H11[H11==Inf] <- .Machine$double.xmax

  H2 <- sapply(u5, function(u){qParetoI(u, b)})
  H22 <- sapply(u6, function(u){qParetoI(u, a)})
  H2[H2==Inf] <- .Machine$double.xmax
  H22[H22==Inf] <- .Machine$double.xmax

  X11 <- R * H1 / H11
  X22 <- R * H2 / H22

  X11[X11==0] <- .Machine$double.xmin
  X22[X22==0] <- .Machine$double.xmin
  X11[X11==Inf] <- .Machine$double.xmax
  X22[X22==Inf] <- .Machine$double.xmax
  
  u <- sapply(X11, function(x){pPPPP(x, al, be, a, b)})
  v <- sapply(X22, function(x){pPPPP(x, al, be, a, b)})
  cbind(u,v)
}


#' Joint CDF of the PPPP model
#'
#' Joint CDF of the PPPP model
#' @param x1,x2:   data input.
#' @param al,be,a,b:   parameters.
#' @keywords Joint CDF
#' @export
#' @examples
#' jpPPPP(2,2, 1.2, 0.8, 1, 1)
jpPPPP <- function(x1,x2,al,be,a,b)
{
  if(x1 <= 0 | x2 <= 0)
  {
    out <- 0
  }else
  {
    xm <- min(x1,x2)
    xp <- max(x1,x2)
    FR <- pPP(xm,al,be)
    tem1 <- a*(xm^(-b)+xp^(-b))/(a+b)
    tem2 <- (a^2)*((xm*xp)^(-b))/((a+b)^2)
    tem3 <- b*(xm^a)/(a+b)
    tem4 <- a*b*(xm^a)*(xp^(-b))/((a+b)^2)
    tem5 <- (b^2)*((xm*xp)^a)/((a+b)^2)
    H1 <- HX_PPPP(b, 0, xm, al, be)
    H2 <- HX_PPPP(2*b, 0, xm, al, be)
    H3 <- HX_PPPP(-a, xm, xp, al, be)
    H4 <- HX_PPPP(b-a, xm, xp, al, be)
    H5 <- HX_PPPP(-2*a, xp, Inf, al, be)
    out <- FR - tem1*H1 + tem2*H2 + tem3*H3 - tem4*H4 + tem5*H5
  }
  return(out)
}

#' Joint CDF of the PPPP copula model
#'
#' Joint CDF of the PPPP copula model
#' @param u,v:   data input.
#' @param al, be, a,b:   parameters.
#' @keywords Joint CDF
#' @export
#' @examples
#' pPPPP_COP(0.9, 0.3, 1.2, 0.8, 1, 1)
pPPPP_COP <- function(u,v,al,be,a,b)
{
  if(u==0 | v==0)
  {
    out <- 0
  }else if(u==1)
  {
    out <- v
  }else if(v==1)
  {
    out <- u 
  }else
  {
    q1 <- qPPPP(u,al,be,a,b)
    q2 <- qPPPP(v,al,be,a,b)
    out <- jpPPPP(q1,q2,al,be,a,b)
  }
  return(out)
}

#' Joint CDF of the PPPP copula model (when a=b=1)
#'
#' Joint CDF of the PPPP copula model (when a=b=1)
#' @param u,v:   data input.
#' @param al, be:   parameters.
#' @keywords Joint CDF
#' @export
#' @examples
#' pPPPP_COP_1(0.9, 0.3, 1.2, 0.8)
pPPPP_COP_1 <- function(u,v,al,be)
{
  pPPPP_COP(u,v,al,be,1,1)
}

#' Joint CDF of the PPPP copula that is rotated for 90 degrees counter clock-wise (when a=b=1)
#'
#' Joint CDF of the PPPP copula that is rotated for 90 degrees counter clock-wise (when a=b=1)
#' @param u,v:   data input.
#' @param al, be:   parameters.
#' @keywords Joint CDF
#' @export
#' @examples
#' pPPPP_COP_1_90(0.9, 0.3, 1.2, 0.8)
pPPPP_COP_1_90 <- function(u,v,al,be)
{
  v - pPPPP_COP(v,1-u,al,be,1,1)
}

intg_spr_PPPP_COP = function(u1u2,al,be,a,b)
{
  u1 <- u1u2[1]
  u2 <- u1u2[2]
  out <- 1 - u1 - u2 + pPPPP_COP(u1,u2,al,be,a,b)
  if(is.finite(out)) out else 0
}

#' Spearman's rho of the PPPP copula
#'
#' Spearman's rho of the PPPP copula
#' @param al,be,a,b:    the two shape parameters.
#' @keywords Spearman's rho
#' @export
#' @examples
#' sprPPPP_COP(3, 6,5.9,6.7)
sprPPPP_COP <- function(al,be,a,b)
{
    tmp <- try(cubature::adaptIntegrate(intg_spr_PPPP_COP, al=al, be=be,a=a,b=b,lowerLimit = c(0,0), upperLimit = c(1,1)), silent = T)
    if(is(tmp,"try-error"))
    {
	    cat("sprPPPP_COP error at: ", al, be, a, b, "\n"); return(NA)
    }else{intg <- tmp$integral}
  out <- 12 * intg - 3
  out
}

intg_tau_PPPP_COP <- function(uv,al,be,a,b)
{
  u <- uv[1]
  v <- uv[2]
  tem1 <- pPPPP_COP(u,v,al,be,a,b)
  tem2 <- dPPPP_COP(u,v,al,be,a,b)
  out <- tem1*tem2
  if(is.finite(out)) out else 0
}

#' Kendall's tau of the PPPP copula
#'
#' Kendall's tau of the PPPP copula
#' @param al,be,a,b    the 4 parameters
#' @keywords Kendall's tau
#' @export
#' @examples
#' tauPPPP_COP(3, 6, 5.9, 6.7)
tauPPPP_COP <- function(al,be,a,b)
{
    tmp <- try(cubature::adaptIntegrate(intg_tau_PPPP_COP, al=al, be=be,a=a,b=b, lowerLimit = c(0,0), upperLimit = c(1,1)), silent = T)
    if(is(tmp,"try-error"))
    {
      cat("tauPPPP_COP error at: ", al, be,a,b, "\n"); return(NA)
    }else{intg <- tmp$integral}
  out <- 4 * intg - 1
  out
}

intg_tau_PPPP_COP_90 <- function(uv,al,be,a,b)
{
  u <- uv[1]
  v <- uv[2]
  tem1 <- v - pPPPP_COP(v,1-u,al,be,a,b)
  tem2 <- dPPPP_COP(v,1-u,al,be,a,b)
  out <- tem1*tem2
  if(is.finite(out)) out else 0
}


#' Kendall's tau of the PPPP copula that is rotated 90 degrees counter clockwise
#'
#' Kendall's tau of the PPPP copula that is rotated 90 degrees counter clockwise
#' @param al,be,a,b    the 4 parameters
#' @keywords Kendall's tau
#' @export
#' @examples
#' tauPPPP_COP_90(3, 6, 5.9, 6.7)
tauPPPP_COP_90 <- function(al,be,a,b,method=1)
{
    tmp <- try(cubature::adaptIntegrate(intg_tau_PPPP_COP_90, al=al, be=be,a=a,b=b, lowerLimit = c(0,0), upperLimit = c(1,1)), silent = T)
    if(is(tmp,"try-error"))
    {
      cat("tauPPPP_COP error at: ", al, be,a,b, "\n"); return(NA)
    }else{intg <- tmp$integral}
  out <- 4 * intg - 1
  out
}
