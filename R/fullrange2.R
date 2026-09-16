#' Uniform score function
#'
#' Transform data to uniform scores, this function is from the R package CopulaModel.
#' @param data Numeric vector, matrix, or data frame. Ties use average ranks. Missing values remain NA and are excluded from the sample size separately for each column.
#' @param aunif Finite plotting-position offset greater than -1. Scores are (rank + aunif)/(n + 1 + 2*aunif), where n counts non-missing values.
#' @returns Uniform scores with the input dimensions and names; data frames return a numeric matrix.
#' @export
#' @keywords distribution
#' @examples
#' data  <- cbind(rnorm(10), rnorm(10))
#' uscore(data)
uscore <- function(data, aunif = -0.5)
{
  if (!is.numeric(aunif) || length(aunif) != 1L ||
      !is.finite(aunif) || aunif <= -1)
    stop("aunif must be a finite scalar greater than -1.")
  score_column <- function(x) {
    if (!is.numeric(x) && !(is.logical(x) && all(is.na(x))))
      stop("data must contain numeric values.")
    out <- rep(NA_real_, length(x))
    names(out) <- names(x)
    observed <- !is.na(x)
    n <- sum(observed)
    if (n > 0L)
      out[observed] <- (rank(x[observed], ties.method = "average") + aunif) /
        (n + 1 + 2 * aunif)
    out
  }
  if (is.null(dim(data))) return(score_column(data))
  if (!is.matrix(data) && !is.data.frame(data))
    stop("data must be a numeric vector, matrix, or data frame.")
  out <- matrix(NA_real_, nrow(data), ncol(data), dimnames = dimnames(data))
  for (j in seq_len(ncol(data))) out[, j] <- score_column(data[, j])
  out
}
#' CDF of univariate margins of the GGEE model
#'
#' CDF of univariate margins of the GGEE model
#' @param x Numeric marginal observations.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param maxit Positive maximum number of numerical iterations. Nonconvergence is reported.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' pGGEE(2, 1, 2)
pGGEE <- function (x, al, be, maxit = 1e+05)
{
    .positive_parameters(al, be)
    if (!is.numeric(x))
        stop("x must be numeric.")
    if (any(!is.finite(x) | x <= 0)) {
        out <- rep(NA_real_, length(x))
        ok <- !is.na(x)
        out[ok & x <= 0] <- 0
        out[ok & x == Inf] <- 1
        inside <- ok & is.finite(x) & x > 0
        out[inside] <- pGGEE(x[inside], al, be, maxit)
        return(out)
    }
    idx1 <- (x >= 100)
    idx2 <- (x <= 100)
    hyp1 <- Re(hypergeo::hypergeo(1, be, al + be + 1, 1 - x[idx1], tol = 1e-12, maxiter = maxit))
    hyp2 <- Re(hypergeo::hypergeo(1, al + 1, al + be + 1, 1 - 1/x[idx2], tol = 1e-12, maxiter = maxit))/x[idx2]
    hyp <- rep(0, length(x))
    hyp[idx1] <- hyp1
    hyp[idx2] <- hyp2
    out <- tryCatch(1 - al/(al + be) * hyp, error = function(err) FALSE, warning = function(err) FALSE)
    if (all(!is.logical(out) & is.finite(out))) {
        return(out)
    }
    else {
        cat("Warning! NA returned!", "\n")
        return(NA)
    }
}
#' Density of univariate margins of the GGEE model
#'
#' Density of univariate margins of the GGEE model
#' @param x Numeric marginal observations.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param maxit Positive maximum number of numerical iterations. Nonconvergence is reported.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' dGGEE(3, 1, 2)
dGGEE <- function (x, al, be, maxit = 1e+05)
{
    .positive_parameters(al, be)
    if (!is.numeric(x))
        stop("x must be numeric.")
    if (any(!is.finite(x) | x <= 0)) {
        out <- rep(NA_real_, length(x))
        ok <- !is.na(x)
        out[ok & (x < 0 | x == Inf)] <- 0
        out[ok & x == 0] <- if (al > 1)
            be/(al - 1)
        else Inf
        inside <- ok & is.finite(x) & x > 0
        out[inside] <- dGGEE(x[inside], al, be, maxit)
        return(out)
    }
    out <- tryCatch((al * be)/(al + be)/(al + be + 1) * Re(hypergeo::hypergeo(2, be + 1, al + be + 2,
        1 - x, tol = 1e-12, maxiter = maxit)), error = function(err) FALSE, warning = function(err) FALSE)
    if (all(!is.logical(out) & is.finite(out)))
        return(out)
    else {
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
#' @param x1 First marginal observation.
#' @param x2 Second marginal observation.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param flag Method flag passed to appell::appellf1().
#' @param integration Logical; use numerical integration instead of Appell's function.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' jdGGEE(10,2, 1, 1)
#' jdGGEE(10,2, 1, 1, integration = TRUE)
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
#' @param u Probability or first copula argument.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param maxit Positive maximum number of numerical iterations. Nonconvergence is reported.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' qGGEE(0.9, 1, 1)
qGGEE <- function(u, al, be, maxit=100000) {
  .positive_parameters(al,be)
  .positive_quantile(u,function(x) pGGEE(x,al,be,maxit),maxit)
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
#' @param n Number of observations to generate.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param seed Optional integer random seed.
#' @param maxit Positive maximum number of numerical iterations. Nonconvergence is reported.
#' @returns A two-column numeric matrix of uniform scores.
#' @export
#' @keywords distribution
#' @examples
#' rGGEE_COP(20, 1.2, 0.2, seed = 100)
rGGEE_COP <- function(n, al, be, seed=NULL, maxit=100000) {
  .positive_parameters(al,be); n<-.sample_size(n)
  if(!is.null(seed)) set.seed(seed)
  lr<-.log_gamma_draw(n,al)-.log_gamma_draw(n,be)
  x<-exp(lr+log(rexp(n))-log(rexp(n)))
  y<-exp(lr+log(rexp(n))-log(rexp(n)))
  cbind(u=pGGEE(x,al,be,maxit),v=pGGEE(y,al,be,maxit))
}
#' Conditional GGEE copula
#'
#' Partial derivative wrt the second argument of the GGEE copula
#' @param u Probability or first copula argument.
#' @param v Conditioning probabilities strictly in (0,1).
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param maxit Positive maximum number of numerical iterations. Nonconvergence is reported.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' C2GGEE_COP(0.2, 0.6, al = 1.2, be = 0.8)
C2GGEE_COP <- function(u,v,al,be,maxit=100000) {
  .positive_parameters(al,be); .probabilities(u); .probabilities(v,TRUE,"v")
  if(length(u)!=1L || length(v)!=1L) stop("u and v must be scalars.")
  if(u==0 || u==1) return(u)
  x1<-qGGEE(u,al,be,maxit); x2<-qGGEE(v,al,be,maxit)
  out<-Dx2_GGEE(x1,x2,al,be)/dGGEE(x2,al,be,maxit)
  if(!is.finite(out) || out<0 || out>1) stop("GGEE conditional CDF evaluation failed.")
  out
}

dGGEE_COP_0 <- function(u,v,al,be,flag=1,integration=FALSE,maxit=100000) {
  .positive_parameters(al,be); .probabilities(c(u,v),TRUE)
  q1<-qGGEE(u,al,be,maxit); q2<-qGGEE(v,al,be,maxit)
  out<-jdGGEE(q1,q2,al,be,flag,integration)/dGGEE(q1,al,be,maxit)/dGGEE(q2,al,be,maxit)
  if(any(!is.finite(out) | out<0)) stop("GGEE density evaluation failed.")
  out
}
#' Copula Density Function - GGEE_COP
#'
#' Copula density function of the bivariate copula that has full-range tail dependence in both upper and lower tails
#' @param u Probability or first copula argument.
#' @param v Second copula argument.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param flag Method flag passed to appell::appellf1().
#' @param integration Logical; use numerical integration instead of Appell's function.
#' @param maxit Positive maximum number of numerical iterations. Nonconvergence is reported.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
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
#' Joint CDF of the GGEE model
#'
#' Joint CDF of the GGEE model
#' @param x1 First marginal observation.
#' @param x2 Second marginal observation.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param flag Method flag passed to appell::appellf1().
#' @param integration Logical; use numerical integration instead of Appell's function.
#' @param maxit Positive maximum number of numerical iterations. Nonconvergence is reported.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' jpGGEE(0.2, 0.4, 1, 1)
#' jpGGEE(0.2, 0.4, 1, 1, integration = TRUE)
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
#' @param u Probability or first copula argument.
#' @param v Second copula argument.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param flag Method flag passed to appell::appellf1().
#' @param integration Logical; use numerical integration instead of Appell's function.
#' @param maxit Positive maximum number of numerical iterations. Nonconvergence is reported.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' pGGEE_COP(0.9, 0.3, 1, 1)
#' pGGEE_COP(0.9, 0.3, 1, 1, integration=TRUE)
pGGEE_COP <- function(u,v,al,be,flag=1,integration=FALSE,maxit=100000) {
  .positive_parameters(al,be); .probabilities(c(u,v))
  if(length(u)!=1L || length(v)!=1L) stop("u and v must be scalars.")
  if(u==0 || v==0) return(0)
  if(u==1) return(v)
  if(v==1) return(u)
  jpGGEE(qGGEE(u,al,be,maxit),qGGEE(v,al,be,maxit),al,be,flag,integration,maxit)
}

intg_tau_E <- function(xy,al,be) {
  # Transform the beta-weighted integral to uniform coordinates.
  # For log odds ratio z, g(z) = exp(-z)*(z-1+exp(-z))/(1-exp(-z))^2.
  # Symmetry gives E[g]=1/2, so tau = E[(2g-1)^2].
  if(any(xy==0 | xy==1)) return(if(xy[1]==xy[2]) 0 else 1)
  x<-qbeta(xy,be,al)
  omx<-qbeta(xy,al,be,lower.tail=FALSE)
  if(any(x<=0 | omx<=0)) stop("GGEE tau exceeds beta-quantile precision; use tauGGEE_COP_sim.")
  z<-abs(diff(log(x)-log(omx)))
  if(z<1e-3) g<-.5-z/6+z^3/180 else {
    q<-exp(-z); g<-q*(z-1+q)/(-expm1(-z))^2
  }
  (2*g-1)^2
}
#' Kendall's tau of the GGEE copula
#'
#' Kendall's tau of the bivariate copula that has full-range tail dependence in both upper and lower tails
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' tauGGEE_COP(0.5, 1)
tauGGEE_COP <- function(al,be) {
  .positive_parameters(al,be)
  .checked_cubature(intg_tau_E,al=al,be=be,lowerLimit=c(0,0),upperLimit=c(1,1))$integral
}
#' Kendall's tau of the GGEE copula (based on simulations)
#'
#' Kendall's tau of the bivariate copula that has full-range tail dependence in both upper and lower tails
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param n Number of observations to generate.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' tauGGEE_COP_sim(0.5, 1, n=10000)
tauGGEE_COP_sim <- function(al,be,n=100000) {
  .positive_parameters(al,be); n<-.sample_size(n)
  if(n<1L) stop("n must be positive.")
  sample_pair<-function() {
    lr<-.log_gamma_draw(n,al)-.log_gamma_draw(n,be)
    cbind(lr+log(rexp(n))-log(rexp(n)),lr+log(rexp(n))-log(rexp(n)))
  }
  x<-sample_pair(); y<-sample_pair()
  mean(sign(x[,1]-y[,1])*sign(x[,2]-y[,2]))
}


intg_spr_C = function(u1u2,al,be,flag=1, integration=F, maxit=100000)
{
  u1 <- u1u2[1]
  u2 <- u1u2[2]
  out <- 1 - u1 - u2 + pGGEE_COP(u1,u2,al,be,flag,integration, maxit=maxit)
  if (!is.finite(out)) stop("Non-finite dependence integrand."); out
}
#' Spearman's rho of the GGEE copula
#'
#' Spearman's rho of the GGEE copula
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param flag Method flag passed to appell::appellf1().
#' @param integration Logical; use numerical integration instead of Appell's function.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' \dontrun{
#' sprGGEE_COP(1.2, 0.6)
#' }
sprGGEE_COP <- function(al,be,flag=1, integration=F)
{
  tmp <- try(.checked_cubature(intg_spr_C, al=al, be=be, flag=flag,integration=integration,lowerLimit = c(0,0), upperLimit = c(1,1)), silent = T)
  if(is(tmp,"try-error"))
  {
   stop(as.character(tmp))
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
#' CDF of univariate margins of the PP model (i.e., Pareto I / Pareto I)
#'
#' CDF of univariate margins of the PP model (i.e., Pareto I / Pareto I)
#' @param x Numeric marginal observations.
#' @param ga1 Positive lower-tail shape parameter.
#' @param ga2 Positive upper-tail shape parameter.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' plot(sapply(seq(0,10, length=50), function(x){pPP(x, 2,2)}), type="l")
pPP <- function(x, ga1, ga2) {
  .positive_parameters(ga1,ga2)
  if (!is.numeric(x)) stop("x must be numeric.")
  out <- rep(NA_real_,length(x)); ok <- !is.na(x)
  out[ok & x<=0] <- 0; out[ok & x==Inf] <- 1
  lo <- ok & x>0 & x<1; hi <- ok & is.finite(x) & x>=1
  out[lo] <- ga2/(ga1+ga2)*x[lo]^ga1
  out[hi] <- -expm1(log(ga1/(ga1+ga2))-ga2*log(x[hi]))
  out
}
#' density of univariate margins of the PP model (i.e., Pareto I / Pareto I)
#'
#' density of univariate margins of the PP model (i.e., Pareto I / Pareto I)
#' @param x Numeric marginal observations.
#' @param ga1 Positive lower-tail shape parameter.
#' @param ga2 Positive upper-tail shape parameter.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' plot(sapply(seq(0,10, length=50), function(x){dPP(x, 2,2)}), type="l")
dPP <- function(x, ga1, ga2) {
  .positive_parameters(ga1,ga2)
  if (!is.numeric(x)) stop("x must be numeric.")
  out <- rep(NA_real_,length(x)); ok <- !is.na(x)
  out[ok & (x<0 | x==Inf)] <- 0
  lo <- ok & x>=0 & x<1; hi <- ok & is.finite(x) & x>=1
  out[lo] <- ga1*ga2/(ga1+ga2)*x[lo]^(ga1-1)
  out[hi] <- ga1*ga2/(ga1+ga2)*x[hi]^(-ga2-1)
  out
}
#' quantile of univariate margins of the PP model (i.e., Pareto I / Pareto I)
#'
#' quantile of univariate margins of the PP model (i.e., Pareto I / Pareto I)
#' @param u Probability or first copula argument.
#' @param ga1 Positive lower-tail shape parameter.
#' @param ga2 Positive upper-tail shape parameter.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' plot(sapply(seq(0.99,1-0.001, length=1000), function(x){qPP(x, 0.2,0.3)}), type="l")
qPP <- function(u, ga1, ga2) {
  .positive_parameters(ga1,ga2); .probabilities(u)
  out <- numeric(length(u)); cutoff <- ga2/(ga1+ga2)
  lower <- u<=cutoff
  out[lower] <- exp((log(u[lower])-log(cutoff))/ga1)
  out[!lower] <- exp((log(ga1/(ga1+ga2))-log1p(-u[!lower]))/ga2)
  out
}

# internal function the g function for pPPPP()
g_PPPP <- function(a1,a2,a3,a4)
{
  1/((1+a1/a2)*(1+a1/a4)*(1-a1/a3))
}
#' CDF of univariate margins of the PPPP model
#'
#' CDF of univariate margins of the PPPP model
#' @param x Numeric marginal observations.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param a Positive shape parameter.
#' @param b Positive shape parameter.
#' @param log Logical; return the logarithm of the CDF.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' plot(sapply(seq(0.0001,20, length=50), function(x){pPPPP(x, 1,1,2,2,log=FALSE)}), type="l",ylab="")
pPPPP <- function(x, al, be, a, b, log=FALSE) {
  .positive_parameters(al,be,a,b)
  if(!is.numeric(x) || length(x)!=1L || is.na(x)) stop("x must be a non-missing numeric scalar.")
  if(!is.logical(log) || length(log)!=1L || is.na(log)) stop("log must be TRUE or FALSE.")
  if(x<=0) return(if(log) -Inf else 0)
  if(x==Inf) return(if(log) 0 else 1)
  FR<-pPP(x,al,be)
  I1<-HX_PPPP(b,0,x,al,be); I2<-HX_PPPP(-a,x,Inf,al,be)
  out<-FR-a/(a+b)*x^(-b)*I1+b/(a+b)*x^a*I2
  if(!is.finite(out) || out<0 || out>1)
    stop("PPPP marginal CDF exceeds numerical precision.")
  if(log) log(out) else out
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
#' @param u Probability or first copula argument.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param a Positive shape parameter.
#' @param b Positive shape parameter.
#' @param maxit Positive maximum number of numerical iterations. Nonconvergence is reported.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' \dontrun{
#' plot(sapply(seq(0.001,0.999, length=100), function(x){qPPPP(x, 0.3, 1.3, 1, 1)}), type="l",ylab="")
#' }
qPPPP <- function(u, al, be, a, b, maxit=100000) {
  .positive_parameters(al,be,a,b)
  .positive_quantile(u,function(x) pPPPP(x,al,be,a,b),maxit)
}
#' Density of univariate margins of the PPPP model
#'
#' Density of univariate margins of the PPPP model
#' @param x Numeric marginal observations.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param a Positive shape parameter.
#' @param b Positive shape parameter.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' plot(sapply(seq(0,20, length=100), function(x){dPPPP(x, 1,1,2,2)}), type="l",ylab="")
dPPPP <- function (x, al, be, a, b)
{
    .positive_parameters(al, be, a, b)
    if (!is.numeric(x) || length(x) != 1L || is.na(x))
        stop("x must be a non-missing numeric scalar.")
    if (x < 0 || x == Inf)
        return(0)
    if (x == 0) {
        m <- min(al, a)
        if (m < 1 || (m == 1 && al == a))
            return(Inf)
        if (m > 1)
            return(0)
        return(if (al == 1) g_PPPP(al, be, a, b) else g_PPPP(a, b, al, be))
    }
    w <- a * b/(a + b)
    HR1 <- HX_PPPP(b, 0, x, al, be)
    HR2 <- HX_PPPP(-a, x, Inf, al, be)
    out <- w * ((x^(-b - 1)) * HR1 + (x^(a - 1)) * HR2)
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
#' @param x1 First marginal observation.
#' @param x2 Second marginal observation.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param a Positive shape parameter.
#' @param b Positive shape parameter.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
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

dPPPP_COP_0 <- function(u,v,al,be,a,b) {
  .positive_parameters(al,be,a,b); .probabilities(c(u,v),TRUE)
  q1<-qPPPP(u,al,be,a,b); q2<-qPPPP(v,al,be,a,b)
  out<-jdPPPP(q1,q2,al,be,a,b)/dPPPP(q1,al,be,a,b)/dPPPP(q2,al,be,a,b)
  if(any(!is.finite(out) | out<0)) stop("PPPP density evaluation failed.")
  out
}
#' Copula Density Function - PPPP_COP
#'
#' PPPP_COP - Copula density function of the bivariate copula that has full-range tail dependence in both upper and lower tails
#' @param u Probability or first copula argument.
#' @param v Second copula argument.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param a Positive shape parameter.
#' @param b Positive shape parameter.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' dPPPP_COP(0.2, 0.4, 1,1,2,2)
dPPPP_COP <- Vectorize(dPPPP_COP_0, c("u", "v"))
#' Copula Density Function of PPPP_COP (when a=b=1)
#'
#' PPPP_COP - Copula density function of PPPP_COP (when a=b=1)
#' @param uvec First copula arguments.
#' @param vvec Second copula arguments.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
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
#' @param uvec First copula arguments.
#' @param vvec Second copula arguments.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
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
#' Conditional PPPP copula
#'
#' Partial derivative wrt the second argument of the PPPP copula
#' @param u Probability or first copula argument.
#' @param v Conditioning probabilities strictly in (0,1).
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param a Positive shape parameter.
#' @param b Positive shape parameter.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' C2PPPP_COP(0.2, 0.6, al=1.2, be=0.8, a=1, b=1)
C2PPPP_COP <- function(u,v,al,be,a,b) {
  .positive_parameters(al,be,a,b); .probabilities(u); .probabilities(v,TRUE,"v")
  if(length(u)!=1L || length(v)!=1L) stop("u and v must be scalars.")
  if(u==0 || u==1) return(u)
  out<-X1onX2_PPPP(qPPPP(u,al,be,a,b),qPPPP(v,al,be,a,b),al,be,a,b)
  if(!is.finite(out) || out<0 || out>1) stop("PPPP conditional CDF evaluation failed.")
  out
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
#' @param n Number of observations to generate.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param a Positive shape parameter.
#' @param b Positive shape parameter.
#' @param seed Optional integer random seed.
#' @returns A two-column numeric matrix of uniform scores.
#' @export
#' @keywords distribution
#' @examples
#' plot(rPPPP_COP(100, 1.2, 0.2, 1, 1, seed = 100))
rPPPP_COP <- function(n,al,be,a,b,seed=NULL) {
  .positive_parameters(al,be,a,b); n<-.sample_size(n)
  if(!is.null(seed)) set.seed(seed)
  u<-matrix(runif(6*n),n,6)
  lr<--log1p(-u[,1])/be+log1p(-u[,2])/al
  lx<-lr-log1p(-u[,3])/b+log1p(-u[,4])/a
  ly<-lr-log1p(-u[,5])/b+log1p(-u[,6])/a
  margin<-function(z) vapply(exp(z),function(x) pPPPP(x,al,be,a,b),numeric(1))
  cbind(u=margin(lx),v=margin(ly))
}
#' Joint CDF of the PPPP model
#'
#' Joint CDF of the PPPP model
#' @param x1 First marginal observation.
#' @param x2 Second marginal observation.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param a Positive shape parameter.
#' @param b Positive shape parameter.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
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
#' @param u Probability or first copula argument.
#' @param v Second copula argument.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param a Positive shape parameter.
#' @param b Positive shape parameter.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' pPPPP_COP(0.9, 0.3, 1.2, 0.8, 1, 1)
pPPPP_COP <- function(u,v,al,be,a,b)
{
  .positive_parameters(al,be,a,b); .probabilities(c(u,v))
  if(length(u)!=1L || length(v)!=1L) stop("u and v must be scalars.")
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
#' @param u Probability or first copula argument.
#' @param v Second copula argument.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' pPPPP_COP_1(0.9, 0.3, 1.2, 0.8)
pPPPP_COP_1 <- function(u,v,al,be)
{
  pPPPP_COP(u,v,al,be,1,1)
}
#' Joint CDF of the PPPP copula that is rotated for 90 degrees counter clock-wise (when a=b=1)
#'
#' Joint CDF of the PPPP copula that is rotated for 90 degrees counter clock-wise (when a=b=1)
#' @param u Probability or first copula argument.
#' @param v Second copula argument.
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
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
  if (!is.finite(out)) stop("Non-finite dependence integrand."); out
}
#' Spearman's rho of the PPPP copula
#'
#' Spearman's rho of the PPPP copula
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param a Positive shape parameter.
#' @param b Positive shape parameter.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' \dontrun{
#' sprPPPP_COP(3, 6,5.9,6.7)
#' }
sprPPPP_COP <- function(al,be,a,b)
{
    tmp <- try(.checked_cubature(intg_spr_PPPP_COP, al=al, be=be,a=a,b=b,lowerLimit = c(0,0), upperLimit = c(1,1)), silent = T)
    if(is(tmp,"try-error"))
    {
	    stop(as.character(tmp))
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
  if (!is.finite(out)) stop("Non-finite dependence integrand."); out
}
#' Kendall's tau of the PPPP copula
#'
#' Kendall's tau of the PPPP copula
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param a Positive shape parameter.
#' @param b Positive shape parameter.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' \dontrun{
#' tauPPPP_COP(3, 6, 5.9, 6.7)
#' }
tauPPPP_COP <- function(al,be,a,b)
{
    tmp <- try(.checked_cubature(intg_tau_PPPP_COP, al=al, be=be,a=a,b=b, lowerLimit = c(0,0), upperLimit = c(1,1)), silent = T)
    if(is(tmp,"try-error"))
    {
      stop(as.character(tmp))
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
  if (!is.finite(out)) stop("Non-finite dependence integrand."); out
}
#' Kendall's tau of the PPPP copula that is rotated 90 degrees counter clockwise
#'
#' Kendall's tau of the PPPP copula that is rotated 90 degrees counter clockwise
#' @param al Positive first shape parameter.
#' @param be Positive second shape parameter.
#' @param a Positive shape parameter.
#' @param b Positive shape parameter.
#' @param method Reserved for compatibility; the current implementation uses numerical integration.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' \dontrun{
#' tauPPPP_COP_90(3, 6, 5.9, 6.7)
#' }
tauPPPP_COP_90 <- function(al,be,a,b,method=1)
{
    tmp <- try(.checked_cubature(intg_tau_PPPP_COP_90, al=al, be=be,a=a,b=b, lowerLimit = c(0,0), upperLimit = c(1,1)), silent = T)
    if(is(tmp,"try-error"))
    {
      stop(as.character(tmp))
    }else{intg <- tmp$integral}
  out <- 4 * intg - 1
  out
}

#########################
# The GGGG copula

## Gauss-Legendre nodes and weights on [0, 1]
gauss_legendre_01 <- function(n = 256)
{
  if (!is.numeric(n) || length(n)!=1L || !is.finite(n) || n<2 || n!=floor(n))
    stop("Quadrature size must be an integer of at least 2.")
  i <- seq_len(n - 1)
  b <- i / sqrt(4 * i^2 - 1)

  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- b
  J[cbind(i + 1, i)] <- b

  eig <- eigen(J, symmetric = TRUE)

  nodes <- eig$values
  weights <- 2 * eig$vectors[1, ]^2

  ord <- order(nodes)

  ## Transform from [-1, 1] to [0, 1]
  q <- (nodes[ord] + 1) / 2
  w <- weights[ord] / 2

  list(nodes = q, weights = w)
}

pGGGG_from_y <- function(x, theta1, theta2, y, one_minus_y, w) {
  .positive_parameters(theta1,theta2)
  if(!is.numeric(x) || anyNA(x)) stop("x must be numeric and non-missing.")
  vapply(x,function(xx) {
    if(xx<=0) return(0)
    if(xx==Inf) return(1)
    z<-plogis(log(xx)+log(one_minus_y)-log(y))
    out<-sum(w*pbeta(z,theta1,theta2))
    if(!is.finite(out)) stop("GGGG marginal CDF exceeds numerical precision.")
    out
  },numeric(1))
}

## Stable marginal CDF pGGGG(x)
pGGGG <- function(x, alpha, beta, theta1, theta2, n_nodes=256, quad=NULL) {
  .positive_parameters(alpha,beta,theta1,theta2)
  if(is.null(quad)) quad<-gauss_legendre_01(n_nodes)
  quad<-.quadrature_rule(quad)
  y<-qbeta(quad$nodes,alpha,beta)
  one_minus_y<-qbeta(quad$nodes,beta,alpha,lower.tail=FALSE)
  pGGGG_from_y(x,theta1,theta2,y,one_minus_y,quad$weights)
}

## ============================================================
## Sampling from the GGGG copula
##
## Model:
## R = R1 / R2
## R1 ~ Gamma(alpha, 1)
## R2 ~ Gamma(beta, 1)
##
## H_i = G_i1 / G_i2
## G_i1 ~ Gamma(theta1, 1)
## G_i2 ~ Gamma(theta2, 1)
##
## X_i = R H_i
## U_i = F_i(X_i)
## ============================================================

rGGGG_COP <- function(n, al, be, th1, th2, seed=NULL, quad=NULL) {
  .positive_parameters(al,be,th1,th2); n<-.sample_size(n)
  if(!is.null(seed)) set.seed(seed)
  if(is.null(quad)) quad<-gauss_legendre_01(512)
  quad<-.quadrature_rule(quad)
  lr<-.log_gamma_draw(n,al)-.log_gamma_draw(n,be)
  lx<-lr+.log_gamma_draw(n,th1)-.log_gamma_draw(n,th2)
  ly<-lr+.log_gamma_draw(n,th1)-.log_gamma_draw(n,th2)
  y<-qbeta(quad$nodes,al,be)
  omy<-qbeta(quad$nodes,be,al,lower.tail=FALSE)
  if(any(y==0 | omy==0)) stop("GGGG beta quadrature exceeds numerical precision for these shapes.")
  margin<-function(z) vapply(z,function(lx)
    sum(quad$weights*pbeta(plogis(lx+log(omy)-log(y)),th1,th2)),numeric(1))
  cbind(u=margin(lx),v=margin(ly))
}


#########################################################################
# FRA1: Full-range tail dependence Archimedean copula 1
# Parameter order: eta (lower tail), theta (upper tail), each in [-1, 1).
#   pFRA1_COP(u, v, eta, theta)       copula CDF
#   dFRA1_COP(u, v, eta, theta)       copula density
#   logdFRA1_COP(u, v, eta, theta)    log copula density
#   C2FRA1_COP(u, v, eta, theta)      dC(u,v)/dv = P(U <= u | V = v)
#   C2invFRA1_COP(p, v, eta, theta)   inverse of C2 with respect to u
#   rFRA1_COP(n, eta, theta, seed)    n-by-2 matrix of uniform scores
#   psiFRA1(t, eta, theta)           generator
#   phiFRA1(u, eta, theta)           inverse generator
#   logderivFRA1(t, eta, theta)      log(-psi') and log(psi'')
#   KFRA1_COP(w, eta, theta)         Kendall distribution
#   tailFRA1_COP(eta, theta)         tail coefficients and orders
#   tauFRA1_COP(eta, theta)          Kendall's tau
#   sprFRA1_COP(eta, theta)          Spearman's rho
#   betaFRA1_COP(eta, theta)         Blomqvist's beta
#   giniFRA1_COP(eta, theta)         Gini's gamma
#   dependenceFRA1_COP(eta, theta)   all four concordance measures


clip01_FRA1 <- function(x, eps = 0) {
  .probabilities(x)
  if (!is.numeric(eps) || length(eps)!=1 || !is.finite(eps) || eps<0 || eps>=.5)
    stop("eps must be a finite scalar in [0,0.5).")
  pmin(pmax(x, eps), 1-eps)
}

recycle_pair_FRA1 <- function(x, y) {
  if (!is.numeric(x) || !is.numeric(y)) stop("Arguments must be numeric.")
  if (!length(x) || !length(y)) return(list(x=numeric(0),y=numeric(0)))
  n <- max(length(x),length(y))
  if (!length(x) %in% c(1,n) || !length(y) %in% c(1,n))
    stop("Arguments must have equal lengths or one must be scalar.")
  list(x=rep_len(x,n),y=rep_len(y,n))
}

logspace_add_FRA1 <- function(log_x, log_y) {
  z <- recycle_pair_FRA1(log_x, log_y)
  log_x <- z$x
  log_y <- z$y
  m <- pmax(log_x, log_y)
  ans <- m + log(exp(log_x - m) + exp(log_y - m))
  both_zero <- is.infinite(log_x) & log_x < 0 &
    is.infinite(log_y) & log_y < 0
  ans[both_zero] <- -Inf
  ans
}

acosh1p_FRA1 <- function(x) {
  # More stable than acosh(1 + x) when x is close to zero.
  log1p(x + sqrt(x) * sqrt(x + 2))
}

coshm1_FRA1 <- function(x) {
  2 * sinh(x / 2)^2
}

check_parameter_FRA1 <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < -1 || x >= 1) {
    stop(name, " must be a finite scalar in [-1, 1).")
  }
  invisible(TRUE)
}

A_FRA1 <- function(s) {
  s <- as.numeric(s)
  if (any(!is.finite(s) | s < 0)) {
    stop("A_fun requires finite nonnegative arguments.")
  }
  acosh1p_FRA1(s)
}

F_FRA1 <- function(s, p) {
  s <- as.numeric(s)
  p <- abs(as.numeric(p)[1])
  if (p > 1) stop("p must be in [0, 1].")
  A <- A_FRA1(s)
  if (p < 1e-7) {
    0.5 * A^2
  } else {
    coshm1_FRA1(p * A) / p^2
  }
}

Fp_FRA1 <- function(s, p) {
  s <- as.numeric(s)
  p <- abs(as.numeric(p)[1])
  if (p > 1) stop("p must be in [0, 1].")

  ans <- numeric(length(s))
  small <- s < 1e-6
  if (any(small)) {
    ss <- s[small]
    ans[small] <- 1 - (1 - p^2) * ss / 3 + (1 - p^2) * (4 - p^2) * ss^2 / 30
  }
  if (any(!small)) {
    ss <- s[!small]
    A <- A_FRA1(ss)
    D <- sqrt(ss) * sqrt(ss + 2)
    if (p < 1e-7) {
      ans[!small] <- A / D
    } else {
      ans[!small] <- sinh(p * A) / (p * D)
    }
  }
  ans
}

Fpp_FRA1 <- function(s, p) {
  s <- as.numeric(s)
  p <- abs(as.numeric(p)[1])
  if (p > 1) stop("p must be in [0, 1].")

  ans <- numeric(length(s))
  small <- s < 1e-5
  if (any(small)) {
    ss <- s[small]
    ans[small] <- -(1 - p^2) / 3 + (1 - p^2) * (4 - p^2) * ss / 15
  }
  if (any(!small)) {
    ss <- s[!small]
    A <- A_FRA1(ss)
    D <- sqrt(ss) * sqrt(ss + 2)
    if (p < 1e-7) {
      ans[!small] <- (D - (1 + ss) * A) / D^3
    } else {
      ans[!small] <-
        (p * cosh(p * A) * D - (1 + ss) * sinh(p * A)) / (p * D^3)
    }
  }
  ans
}

Finv_FRA1 <- function(x, p) {
  x <- as.numeric(x)
  p <- abs(as.numeric(p)[1])
  if (any(!is.finite(x) | x < 0)) {
    stop("Finv_fun requires finite nonnegative arguments.")
  }
  if (p > 1) stop("p must be in [0, 1].")

  if (p < 1e-7) {
    coshm1_FRA1(sqrt(2 * x))
  } else {
    coshm1_FRA1(acosh1p_FRA1(p^2 * x) / p)
  }
}

d_FRA1 <- function(p) {
  Finv_FRA1(1, p)
}

S_bundle_FRA1 <- function(x, p) {
  x <- as.numeric(x)
  if (any(!is.finite(x) | x < 0)) {
    stop("S_bundle requires finite nonnegative arguments.")
  }
  p <- abs(as.numeric(p)[1])
  d <- d_FRA1(p)

  T <- Finv_FRA1(F_FRA1(x, p) + 1, p)
  S <- T - d

  # T-d loses precision when x is very small. Use the local Taylor expansion
  # S(x) = T'(0)x + T''(0)x^2/2 + o(x^2) in that region.
  small <- x < 1e-5 * (1 + d)
  fp_d <- Fp_FRA1(d, p)
  fpp_d <- Fpp_FRA1(d, p)
  t1_0 <- 1 / fp_d
  t2_0 <- Fpp_FRA1(0, p) / fp_d - fpp_d / fp_d^3
  if (any(small)) {
    S[small] <- t1_0 * x[small] + 0.5 * t2_0 * x[small]^2
    T[small] <- d + S[small]
  }

  fp_x <- Fp_FRA1(x, p)
  fp_T <- Fp_FRA1(T, p)
  fpp_x <- Fpp_FRA1(x, p)
  fpp_T <- Fpp_FRA1(T, p)

  Sp <- fp_x / fp_T
  Spp <- fpp_x / fp_T - fp_x^2 * fpp_T / fp_T^3

  list(S = S, Sp = Sp, Spp = Spp, T = T, d = d)
}

J_bundle_FRA1 <- function(t, p) {
  t <- as.numeric(t)
  if (any(!is.finite(t) | t <= 0)) {
    stop("J_bundle requires finite positive arguments.")
  }
  x <- 1 / t
  sb <- S_bundle_FRA1(x, p)

  S <- sb$S
  Sp <- sb$Sp
  Spp <- sb$Spp

  J <- 1 / S
  Jp <- Sp / (t^2 * S^2)

  term1 <- -2 * Sp / (t^3 * S^2)
  term2 <- -Spp / (t^4 * S^2)
  term3 <- 2 * Sp^2 / (t^4 * S^3)
  Jpp <- term1 + term2 + term3

  # Remove only tiny positive roundoff errors. A materially positive value is
  # retained and will be rejected by the density evaluator.
  scale <- abs(term1) + abs(term2) + abs(term3) + .Machine$double.xmin
  roundoff <- Jpp > 0 & Jpp <= 1e-8 * scale
  Jpp[roundoff] <- 0

  list(value = J, first = Jp, second = Jpp)
}

Jinv_FRA1 <- function(y, p) {
  y <- as.numeric(y)
  if (any(!is.finite(y) | y <= 0)) {
    stop("Jinv_fun requires finite positive arguments.")
  }
  p <- abs(as.numeric(p)[1])
  d <- d_FRA1(p)
  z <- 1 / y

  delta_F <- F_FRA1(d + z, p) - 1
  small <- z < 1e-5 * (1 + d)
  if (any(small)) {
    delta_F[small] <- Fp_FRA1(d, p) * z[small] + 0.5 * Fpp_FRA1(d, p) * z[small]^2
  }
  1 / Finv_FRA1(delta_F, p)
}

K_bundle_FRA1 <- function(t, a) {
  t <- as.numeric(t)
  a <- as.numeric(a)[1]
  if (a <= 0 || a > 1) stop("a must be in (0, 1].")
  if (any(!is.finite(t) | t <= 0)) {
    stop("K_bundle requires finite positive arguments.")
  }

  ta <- t^a
  value <- expm1(log1p(ta) / a)
  first <- t^(a - 1) * (1 + ta)^(1 / a - 1)
  second <- (a - 1) * t^(a - 2) * (1 + ta)^(1 / a - 2)
  list(value = value, first = first, second = second)
}

Kinv_FRA1 <- function(s, a) {
  s <- as.numeric(s)
  a <- as.numeric(a)[1]
  if (a <= 0 || a > 1) stop("a must be in (0, 1].")
  if (any(!is.finite(s) | s < 0)) {
    stop("Kinv_fun requires finite nonnegative arguments.")
  }
  exp(log(expm1(a * log1p(s))) / a)
}

W_bundle_FRA1 <- function(t, theta) {
  check_parameter_FRA1(theta, "theta")
  t <- as.numeric(t)
  if (any(!is.finite(t) | t <= 0)) {
    stop("W_bundle requires finite positive arguments.")
  }

  if (theta <= 0) {
    J_bundle_FRA1(t, -theta)
  } else {
    a <- 1 - theta
    kb <- K_bundle_FRA1(t, a)
    jb <- J_bundle_FRA1(kb$value, 0)
    list(
      value = jb$value,
      first = jb$first * kb$first,
      second = jb$second * kb$first^2 + jb$first * kb$second
    )
  }
}

Winv_FRA1 <- function(y, theta) {
  check_parameter_FRA1(theta, "theta")
  y <- as.numeric(y)
  if (theta <= 0) {
    Jinv_FRA1(y, -theta)
  } else {
    Kinv_FRA1(Jinv_FRA1(y, 0), 1 - theta)
  }
}

b_FRA1 <- function(eta) {
  eta <- as.numeric(eta)[1]
  if (eta <= 0 || eta >= 1) stop("b_eta requires eta in (0, 1).")
  (1 - eta) / eta^2
}

Ginv_FRA1 <- function(u, eta, eps = 1e-12) {
  check_parameter_FRA1(eta, "eta")
  u <- clip01_FRA1(as.numeric(u), eps)
  if (eta <= 0) {
    Finv_FRA1(-log(u), -eta)
  } else {
    b <- b_FRA1(eta)
    exponent <- -log(u) / b
    if (any(exponent > 700)) {
      stop("G inverse overflow: move eta farther from 1 or clip the data more.")
    }
    x <- b * expm1(exponent)
    Finv_FRA1(x, eta)
  }
}

G_terms_FRA1 <- function(s, eta) {
  check_parameter_FRA1(eta, "eta")
  s <- as.numeric(s)
  if (any(!is.finite(s) | s < 0)) {
    stop("outer_terms requires finite nonnegative arguments.")
  }

  p <- abs(eta)
  f0 <- F_FRA1(s, p)
  f1 <- Fp_FRA1(s, p)
  f2 <- Fpp_FRA1(s, p)

  if (eta <= 0) {
    log_G <- -f0
    rate <- f1                         # -G'/G
    curvature <- f1^2 - f2             # G''/G
  } else {
    b <- b_FRA1(eta)
    log_Q <- log1p(f0 / b)
    inv_Q <- exp(-log_Q)
    log_G <- -b * log_Q
    rate <- f1 * inv_Q
    curvature <- -f2 * inv_Q + ((b + 1) / b) * (f1 * inv_Q)^2
  }

  list(log_G = log_G, rate = rate, curvature = curvature)
}
#' Generator of the FRA1 copula
#'
#' @param t Nonnegative generator argument; derivatives require finite positive values.
#' @param eta Lower-tail parameter in [-1,1).
#' @param theta Upper-tail parameter in [-1,1).
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
psiFRA1 <- function(t, eta, theta) {
  check_parameter_FRA1(eta,"eta"); check_parameter_FRA1(theta,"theta")
  if (!is.numeric(t) || any(is.na(t) | t<0)) stop("t must be nonnegative.")
  ans <- numeric(length(t)); ans[t==0] <- 1
  inside <- is.finite(t) & t>0
  ans[inside] <- exp(.fra_logderiv(log(t[inside]),eta,theta)$value)
  ans
}
#' Inverse generator of the FRA1 copula
#'
#' @param u Probability or first copula argument.
#' @param eta Lower-tail parameter in [-1,1).
#' @param theta Upper-tail parameter in [-1,1).
#' @param eps Optional clipping threshold in [0,0.5). Default 0 preserves the input; positive values explicitly clip probabilities to [eps,1-eps].
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
phiFRA1 <- function(u, eta, theta, eps = 0) {
  check_parameter_FRA1(eta,"eta"); check_parameter_FRA1(theta,"theta")
  u <- clip01_FRA1(u,eps)
  ans <- numeric(length(u)); ans[u==0] <- Inf
  inside <- u>0 & u<1
  ans[inside] <- exp(.fra_logphi(u[inside],eta,theta))
  if (any(!is.finite(ans[inside])))
    stop("Inverse generator exceeds double precision; use copula functions directly.")
  ans
}
#' Log derivatives of the FRA1 generator
#'
#' @param t Nonnegative generator argument; derivatives require finite positive values.
#' @param eta Lower-tail parameter in [-1,1).
#' @param theta Upper-tail parameter in [-1,1).
#' @returns A list with value (log generator), log_negative_first, and log_second.
#' @export
#' @keywords distribution
logderivFRA1 <- function(t, eta, theta) {
  check_parameter_FRA1(eta,"eta"); check_parameter_FRA1(theta,"theta")
  if (!is.numeric(t) || any(!is.finite(t) | t<=0)) stop("t must be finite and positive.")
  .fra_logderiv(log(t),eta,theta)
}
#' Log-density of the FRA1 copula
#'
#' @param u First probabilities, strictly in (0,1) unless positive eps is supplied.
#' @param v Second probabilities, strictly in (0,1) unless positive eps is supplied.
#' @param eta Lower-tail parameter in [-1,1).
#' @param theta Upper-tail parameter in [-1,1).
#' @param eps Optional clipping threshold in [0,0.5). Default 0 preserves the input; positive values explicitly clip probabilities to [eps,1-eps].
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' logdFRA1_COP(0.3, 0.4, eta=0.2, theta=0.3)
logdFRA1_COP <- function(u, v, eta, theta, eps = 0) {
  check_parameter_FRA1(eta,"eta"); check_parameter_FRA1(theta,"theta")
  z <- recycle_pair_FRA1(u,v)
  u <- clip01_FRA1(z$x,eps); v <- clip01_FRA1(z$y,eps)
  .probabilities(u,TRUE); .probabilities(v,TRUE)
  x <- .fra_logphi(u,eta,theta); y <- .fra_logphi(v,eta,theta)
  ans <- .fra_logderiv(logspace_add_FRA1(x,y),eta,theta)$log_second -
    .fra_logderiv(x,eta,theta)$log_negative_first -
    .fra_logderiv(y,eta,theta)$log_negative_first
  if (any(!is.finite(ans))) stop("FRA1 log-density exceeds numerical precision.")
  ans
}
#' Density of the FRA1 copula
#'
#' Density of the FRA1 copula
#' @param u First probabilities, strictly in (0,1) unless positive eps is supplied.
#' @param v Second probabilities, strictly in (0,1) unless positive eps is supplied.
#' @param eta Lower-tail parameter in [-1,1).
#' @param theta Upper-tail parameter in [-1,1).
#' @param eps Optional clipping threshold in [0,0.5). Default 0 preserves the input; positive values explicitly clip probabilities to [eps,1-eps].
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' dFRA1_COP(c(0.2, 0.3), c(0.5, 0.8), eta = -0.2, theta = 0.3)
#'
dFRA1_COP <- function(u, v, eta, theta, eps = 0) {
  exp(logdFRA1_COP(u,v,eta,theta,eps))
}
#' CDF of the FRA1 copula
#'
#' @param u First probabilities in [0,1]; endpoints are exact when eps is zero.
#' @param v Second probabilities in [0,1]; endpoints are exact when eps is zero.
#' @param eta Lower-tail parameter in [-1,1).
#' @param theta Upper-tail parameter in [-1,1).
#' @param eps Optional clipping threshold in [0,0.5). Default 0 preserves the input; positive values explicitly clip probabilities to [eps,1-eps].
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' pFRA1_COP(0.3, 0.4, eta=0.2, theta=0.3)
pFRA1_COP <- function(u, v, eta, theta, eps = 0) {
  check_parameter_FRA1(eta,"eta"); check_parameter_FRA1(theta,"theta")
  z <- recycle_pair_FRA1(u,v)
  u <- clip01_FRA1(z$x,eps); v <- clip01_FRA1(z$y,eps)
  ans <- pmin(u,v)
  inside <- u>0 & u<1 & v>0 & v<1
  x <- .fra_logphi(u[inside],eta,theta); y <- .fra_logphi(v[inside],eta,theta)
  ans[inside] <- exp(.fra_logderiv(logspace_add_FRA1(x,y),eta,theta)$value)
  if (any(!is.finite(ans))) stop("FRA1 CDF exceeds numerical precision.")
  ans
}
#' Conditional CDF of the FRA1 copula
#'
#' @param u Probability or first copula argument.
#' @param v Conditioning probabilities strictly in (0,1).
#' @param eta Lower-tail parameter in [-1,1).
#' @param theta Upper-tail parameter in [-1,1).
#' @param eps Optional clipping threshold in [0,0.5). Default 0 preserves the input; positive values explicitly clip probabilities to [eps,1-eps].
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
C2FRA1_COP <- function(u, v, eta, theta, eps = 0) {
  check_parameter_FRA1(eta,"eta"); check_parameter_FRA1(theta,"theta")
  z <- recycle_pair_FRA1(u,v)
  u <- clip01_FRA1(z$x,eps); v <- clip01_FRA1(z$y,eps)
  .probabilities(v,TRUE,"v")
  ans <- u; inside <- u>0 & u<1
  x <- .fra_logphi(u[inside],eta,theta); y <- .fra_logphi(v[inside],eta,theta)
  ans[inside] <- exp(.fra_logderiv(logspace_add_FRA1(x,y),eta,theta)$log_negative_first -
    .fra_logderiv(y,eta,theta)$log_negative_first)
  if (any(!is.finite(ans))) stop("FRA1 conditional CDF exceeds numerical precision.")
  pmin(pmax(ans,0),1)
}
#' Conditional quantiles of the FRA1 copula
#'
#' @param p Conditional probabilities in [0,1].
#' @param v Conditioning probabilities strictly in (0,1).
#' @param eta Lower-tail parameter in [-1,1).
#' @param theta Upper-tail parameter in [-1,1).
#' @param eps Optional clipping threshold in [0,0.5). Default 0 preserves the input; positive values explicitly clip probabilities to [eps,1-eps].
#' @param root_tol Positive tolerance for inversion on the log-generator scale.
#' @param max_upper Optional positive upper bound on the generator argument; Inf imposes no finite bound.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
C2invFRA1_COP <- function(p, v, eta, theta, eps=0, root_tol=1e-10, max_upper=Inf) {
  check_parameter_FRA1(eta,"eta"); check_parameter_FRA1(theta,"theta")
  if(length(root_tol)!=1L || !is.finite(root_tol) || root_tol<=0)
    stop("root_tol must be a finite positive scalar.")
  if(length(max_upper)!=1L || is.na(max_upper) || max_upper<=0)
    stop("max_upper must be a positive scalar, possibly Inf.")
  z<-recycle_pair_FRA1(p,v)
  p<-clip01_FRA1(z$x,eps); v<-clip01_FRA1(z$y,eps)
  .probabilities(v,TRUE,"v")
  ans<-p; inside<-which(p>0 & p<1)
  for(i in inside) {
    ly<-.fra_logphi(v[i],eta,theta)
    ld<-.fra_logderiv(ly,eta,theta)$log_negative_first
    residual<-function(lx) {
      out<-.fra_logderiv(logspace_add_FRA1(lx,ly),eta,theta)$log_negative_first-ld-log(p[i])
      if(!is.finite(out)) stop("Inverse conditional CDF exceeds numerical precision.")
      out
    }
    lo<-hi<-min(ly,log(max_upper)); step<-2
    while(residual(lo)<0) {lo<-lo-step; step<-2*step; if(!is.finite(lo)) stop("Could not bracket conditional quantile.")}
    step<-2
    while(residual(hi)>0 && hi<log(max_upper)) {
      hi<-min(hi+step,log(max_upper)); step<-2*step
      if(!is.finite(hi)) stop("Could not bracket conditional quantile.")
    }
    if(residual(hi)>0) stop("Conditional quantile exceeds max_upper.")
    root<-uniroot(residual,c(lo,hi),tol=root_tol,maxiter=1000,check.conv=TRUE)$root
    ans[i]<-exp(.fra_logderiv(root,eta,theta)$value)
  }
  actual<-C2FRA1_COP(ans,v,eta,theta)
  if(any(!is.finite(actual) | abs(actual-p)>max(1e-7,10*root_tol)))
    stop("Conditional quantile cannot be represented at the requested precision.")
  ans
}
#' Kendall distribution of the FRA1 copula
#'
#' @param w Probabilities in [0,1].
#' @param eta Lower-tail parameter in [-1,1).
#' @param theta Upper-tail parameter in [-1,1).
#' @param eps Optional clipping threshold in [0,0.5). Default 0 preserves the input; positive values explicitly clip probabilities to [eps,1-eps].
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
KFRA1_COP <- function(w, eta, theta, eps = 0) {
  check_parameter_FRA1(eta,"eta"); check_parameter_FRA1(theta,"theta")
  w <- clip01_FRA1(w,eps)
  ans <- w; inside <- w>0 & w<1
  lp <- .fra_logphi(w[inside],eta,theta)
  ans[inside] <- w[inside] + exp(lp+.fra_logderiv(lp,eta,theta)$log_negative_first)
  pmin(pmax(ans,0),1)
}
#' Tail coefficients and orders of the FRA1 copula
#'
#' @param eta Lower-tail parameter in [-1,1).
#' @param theta Upper-tail parameter in [-1,1).
#' @returns A named numeric vector of dependence summaries.
#' @export
#' @keywords distribution
tailFRA1_COP <- function(eta, theta) {
  check_parameter_FRA1(theta, "theta")
  check_parameter_FRA1(eta, "eta")

  if (theta <= 0) {
    lambda_U <- 0
    kappa_U <- 1 - theta
  } else {
    lambda_U <- 2 - 2^(1 - theta)
    kappa_U <- 1
  }

  if (eta <= 0) {
    lambda_L <- 0
    kappa_L <- 2^(-eta)
  } else {
    lambda_L <- 2^(-(1 - eta) / eta)
    kappa_L <- 1
  }

  c(lambda_U = lambda_U, kappa_U = kappa_U,
    lambda_L = lambda_L, kappa_L = kappa_L)
}
#' Generate random samples of the FRA1 copula
#'
#' Generate random samples of the FRA1 copula
#' @param n Number of observations to generate.
#' @param eta Lower-tail parameter in [-1,1).
#' @param theta Upper-tail parameter in [-1,1).
#' @param seed Optional integer random seed.
#' @param root_tol Positive tolerance for inversion on the log-generator scale.
#' @param max_upper Optional positive upper bound on the generator argument; Inf imposes no finite bound.
#' @returns A two-column numeric matrix of uniform scores.
#' @export
#' @keywords distribution
#' @examples
#' rFRA1_COP(100, eta = -0.2, theta = 0.3, seed = 123)
#'
rFRA1_COP <- function(n, eta, theta, seed=NULL, root_tol=1e-9, max_upper=Inf) {
  check_parameter_FRA1(eta,"eta"); check_parameter_FRA1(theta,"theta")
  n<-.sample_size(n)
  if(!is.null(seed)) set.seed(seed)
  u<-pmax(runif(n),.Machine$double.xmin); r<-runif(n)
  v<-C2invFRA1_COP(r,u,eta,theta,root_tol=root_tol,max_upper=max_upper)
  cbind(u=u,v=v)
}

clean_concordance_FRA1 <- function(x, name, tolerance = 5e-6) {
  if (!is.finite(x)) return(NA_real_)
  if (x < -tolerance || x > 1 + tolerance) {
    warning(name, " = ", signif(x, 7),
            " lies outside [0,1]; increase quad_n or inspect the parameters.")
    return(x)
  }
  pmin(pmax(x, 0), 1)
}
#' Blomqvist's beta of the FRA1 copula
#'
#' Blomqvist's beta of the FRA1 copula
#' @param eta Lower-tail parameter in [-1,1).
#' @param theta Upper-tail parameter in [-1,1).
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' betaFRA1_COP(eta = -0.2, theta = 0.3)
#'
betaFRA1_COP <- function(eta, theta) {
  4 * pFRA1_COP(0.5, 0.5, eta, theta) - 1
}
#' Dependence measures of the FRA1 copula
#'
#' Calculate Kendall's tau, Spearman's rho, Blomqvist's beta,
#' and Gini's gamma.
#'
#' @param eta Lower-tail parameter in [-1,1).
#' @param theta Upper-tail parameter in [-1,1).
#' @param quad_n Number of Gauss-Legendre nodes, an integer of at least 2.
#' @param quadrature_rule Optional list of interior nodes and positive weights summing to one.
#' @returns A named numeric vector of dependence summaries.
#' @export
#' @keywords distribution
#' @examples
#' dependenceFRA1_COP(eta = -0.2, theta = 0.3)
dependenceFRA1_COP <- function(eta, theta, quad_n=48, quadrature_rule=NULL) {
  check_parameter_FRA1(eta,"eta"); check_parameter_FRA1(theta,"theta")
  if(is.null(quadrature_rule)) quadrature_rule<-gauss_legendre_01(quad_n)
  quadrature_rule<-.quadrature_rule(quadrature_rule)
  u<-quadrature_rule$nodes; weight<-quadrature_rule$weights
  lp<-.fra_logphi(u,eta,theta)
  ld<-.fra_logderiv(lp,eta,theta)$log_negative_first
  tau<-1-4*sum(weight*exp(lp+ld))
  sum_phi<-outer(lp,lp,function(x,y) logspace_add_FRA1(x,y))
  C<-matrix(exp(.fra_logderiv(as.vector(sum_phi),eta,theta)$value),length(u),length(u))
  rho<-12*sum((weight %o% weight)*C)-3
  beta<-betaFRA1_COP(eta,theta)
  gamma<-4*sum(weight*(pFRA1_COP(u,u,eta,theta)+pFRA1_COP(u,1-u,eta,theta)))-2
  c(kendall_tau=clean_concordance_FRA1(tau,"Kendall's tau"),
    spearman_rho=clean_concordance_FRA1(rho,"Spearman's rho"),
    blomqvist_beta=clean_concordance_FRA1(beta,"Blomqvist's beta"),
    gini_gamma=clean_concordance_FRA1(gamma,"Gini's gamma"))
}
#' Kendall's tau of the FRA1 copula
#'
#' Kendall's tau of the FRA1 copula
#' @param eta Lower-tail parameter in [-1,1).
#' @param theta Upper-tail parameter in [-1,1).
#' @param quad_n Number of Gauss-Legendre nodes, an integer of at least 2.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' tauFRA1_COP(eta = -0.2, theta = 0.3)
#'
tauFRA1_COP <- function(eta, theta, quad_n = 64) {
  unname(dependenceFRA1_COP(eta, theta, quad_n)["kendall_tau"])
}
#' Spearman's rho of the FRA1 copula
#'
#' Spearman's rho of the FRA1 copula
#' @param eta Lower-tail parameter in [-1,1).
#' @param theta Upper-tail parameter in [-1,1).
#' @param quad_n Number of Gauss-Legendre nodes, an integer of at least 2.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' sprFRA1_COP(eta = -0.2, theta = 0.3)
#'
sprFRA1_COP <- function(eta, theta, quad_n = 64) {
  unname(dependenceFRA1_COP(eta, theta, quad_n)["spearman_rho"])
}
#' Gini's gamma of the FRA1 copula
#'
#' Gini's gamma of the FRA1 copula
#' @param eta Lower-tail parameter in [-1,1).
#' @param theta Upper-tail parameter in [-1,1).
#' @param quad_n Number of Gauss-Legendre nodes, an integer of at least 2.
#' @returns Numeric function values. Invalid parameters or numerical failures are reported.
#' @export
#' @keywords distribution
#' @examples
#' giniFRA1_COP(eta = -0.2, theta = 0.3)
#'
giniFRA1_COP <- function(eta, theta, quad_n = 64) {
  unname(dependenceFRA1_COP(eta, theta, quad_n)["gini_gamma"])
}
