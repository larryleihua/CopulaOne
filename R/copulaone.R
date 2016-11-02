# copulaone.R

# used for parallel computing
seqRun <- function(i, dat, nco, par, flag=1)
{
  # dataset is divided in nco parts, and the ith part is working
  dat <- as.matrix(dat)
  N <- dim(dat)[1]
  Ni <- floor( N / nco )
  if(i < nco)
  {
    dat_i <- dat[ ((i-1)*Ni+1) : (i*Ni) , ]
  }else if(i == nco) # the last one has all the rest data, maybe greater than Ni
  {
    dat_i <- dat[ ((nco-1)*Ni+1) : N , ]
  }
  den <- dGGEE_COP(u=dat_i[,1], v=dat_i[,2], a=par[1], b=par[2], flag=flag)
  return(-log(den))
}


#' Model fitting with CopulaOne
#'
#' Fit bivariate data with the full-range tail dependence copula based on maximum likelihood
#' @param par: initial parameters for (a, b)
#' @param dat: input of data.
#' @param flag: indicate which numerical method for the appell function (default: flag = 1)
#' @param se: whether standard errors of parameters are reported (default: se = F)
#' @param trace, printlevel: integers for optim() and nlm(), respectively. (default: 0,0)
#' @keywords fitting
#' @export
#' @examples
#' # DO NOT run, and it takes time!
#' data("euro0306")
#' dat <- uscore(euro0306[,c(2,3)])[1:100,]
#' par <- c(0.3, 0.3)
#' fit <- fitCopulaOne(par, dat)
fitCopulaOne <- function(par, dat, flag=1, opt="L-BFGS-B", se=F, lower=c(0.1, 0.1), upper=c(5, 5), trace=0, factr=1e9, printlevel=0)
{
  dat <- as.matrix(dat)
  
  # check positive / negative dependence by 
  ken <- cor(x=dat[,1], y=dat[,2], method = "kendall")
  if(ken < 0){dat[,2] <- 1 - dat[,2]}
  
  nco <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(nco)
  parallel::clusterEvalQ(cl, library(CopulaOne))
  
  if(opt == "nlm")
  {
    obj <- function(par)
    {
      par <- exp(par)+0.01 # parameters are both > 0
      nllk_each_i <- try(parallel::parLapply(cl, seq_len(nco), seqRun, dat=dat, nco=nco, par=par, flag=flag), silent = T) 
      if(is(nllk_each_i,"try-error")){return(20)}else{
        nllk_each_i <- unlist(nllk_each_i)
        out <- sum(nllk_each_i[is.finite(nllk_each_i)])
        return(out)
      }
    }
    if(se==T){hes <- T}else{hes <- F}
    fit <- nlm(obj, p=log(par-0.01), print.level=printlevel, hessian = hes)
    parallel::stopCluster(cl)
    fit$estimate <- exp(fit$estimate)+0.01
    if(se==T)
    {
      va <- diag(solve(fit$hessian))
      err <- rep(NA, length(par)) 
      err[va>0] <- sqrt(va[va>0])*(fit$estimate[va>0]-0.01)
      fit$err <- err
    }
  }else if(opt == "L-BFGS-B")
  {
    obj <- function(par)
    {
      nllk_each_i <- try(parallel::parLapply(cl, seq_len(nco), seqRun, dat=dat, nco=nco, par=par, flag=flag), silent = T) 
      if(is(nllk_each_i,"try-error")){return(20)}else{
        nllk_each_i <- unlist(nllk_each_i)
        out <- sum(nllk_each_i[is.finite(nllk_each_i)])
        return(out)
      }
    }
    if(se==T){hes <- T}else{hes <- F}
    fit <- optim(par=par, obj, method="L-BFGS-B", control=list(trace=trace, factr=factr), hessian=hes, lower=lower, upper=upper)
    parallel::stopCluster(cl)
    if(se==T)
    {
      va <- diag(solve(fit$hessian))
      err <- rep(NA, length(par)) 
      err[va>0] <- sqrt(va[va>0])
      fit$err <- err
    }
  }
  
  if(ken < 0)
  {
    cat("The second variable was transformed to 1 - v", "\n")
    fit$warning = "The second variable was transformed to 1 - v"
  }
  fit
}


#' Contour plots of CopulaOne
#'
#' Contour plots based on either normal scores or uniform scores
#' @param a, b: parameters
#' @param marg: indicate normal scores or uniform scores (default: marg = "normal")
#' @param resolution: number of evaluations for each dimension (default is 30)
#' @param flag: used in the appell function (default is 1)
#' @keywords contour plot
#' @export
#' @examples 
#' plotCopulaOne(a=0.5, b=1.8)
#' plotCopulaOne(a=0.5, b=1.8, marg="uniform", resolution=20)
plotCopulaOne <- function(a, b, marg="normal", flag=1, resolution=30)
{
  zvec <- seq(-2.5, 2.5, length=resolution)
  f <- dnorm(zvec)
  Fvec <- pnorm(zvec)
  nn <- length(zvec)

  nco <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(nco)
  parallel::clusterEvalQ(cl, library(CopulaOne))
  
  Fmat <- cbind(rep(Fvec, each = nn), rep(Fvec, times = nn))
  
  denvec <- try(parallel::parLapply(cl, seq_len(nco), seqRun, dat=Fmat, nco=nco, par=c(a,b), flag=flag), silent = T)
  parallel::stopCluster(cl)
  
  if(is(denvec,"try-error")){cat("density calculated error!", "\n"); return(NA)}else{
    denvec <- exp(-unlist(denvec))
    denvec[!is.finite(denvec)] <- 0
  }
  
  if(marg == "normal")
  {
    denvec <- denvec * rep(f, each = nn) * rep(f, times = nn)
    denmat <- matrix(denvec, nn, nn)
    rangeforplot <- denvec[denvec>max(denmat[1,], denmat[,1], denmat[nn,], denmat[,nn])]
    contour(zvec,zvec, denmat, drawlabels=F, levels = quantile(rangeforplot,seq(0.05, 0.95, length=10), na.rm=T), labcex = 1, main=paste("Contour plot of normal scores", " (a=", format(a), ", b=", format(b), ")", sep=""), xlim = c(-2.5,2.5), ylim = c(-2.5,2.5))
  } else if(marg == "uniform")
  {
    denmat <- matrix(denvec, nn, nn)
    rangeforplot <- denvec
    contour(Fvec, Fvec, denmat, drawlabels=F, levels = quantile(rangeforplot,seq(0.05, 0.95, length=10), na.rm=T), main=paste("Contour plot of uniform scores", " (a=", format(a), ", b=", format(b), ")", sep=""), xlim = c(0,1), ylim = c(0,1))
  }
  invisible(denvec)
}
