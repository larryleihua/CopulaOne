# copulaone.R

# used for parallel computing
seqRun <- function(i, dat, nco, para, flag=1, integration=F, copula_family="PPPP")
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
  if(copula_family=="PPPP")
  {
    den <- dPPPP_COP(u=dat_i[,1], v=dat_i[,2], al=para[1], be=para[2], a=para[3], b=para[4])
  }else if(copula_family=="GGEE")
  {
    den <- dGGEE_COP(u=dat_i[,1], v=dat_i[,2], a=para[1], b=para[2], flag=flag, integration=integration)
  }else
  {
    cat("Please specify copula families from (GGEE, PPPP)", "\n")
    return(1)
  }
  return(-log(den))
}

#' Model fitting with CopulaOne
#'
#' Fit bivariate data with the full-range tail dependence copula based on maximum likelihood
#' @param par: initial parameters for (a, b)
#' @param patternpar: a vector for assigning parameters to be estimated, and this allows same values of different parameters. For example, if patternpar=c(1,2,3,0,3), then there are 3 parameters to estimate, the 3rd and 5th ones are the same, and the 4th parameter is fixed as specified in the initial parameters (i.e., par)
#' @param dat: input of data.
#' @param flag: indicate which numerical method for the appell function (default: flag = 1)
#' @param integration: (Experimental!) using integration instead of appellf1()
#' @param se: whether standard errors of parameters are reported (default: se = F)
#' @param trace, printlevel: integers for optim() and nlm(), respectively. (default: 0,0)
#' @keywords fitting
#' @export
#' @examples
#' \dontrun{
#' data("euro0306")
#' dat <- uscore(euro0306[,c(2,3)])[1:100,]
#' par0 <- c(0.3, 0.3)
#' fit <- fitCopulaOne(par0, dat=dat, copula_family="GGEE")
#' par0 <- c(0.3, 0.3, 1, 1)
#' lower <- rep(0.1, 2)
#' upper <- rep(5, 2)
#' par0 <- c(0.3, 0.3, 1, 1)
#' patternpar <- c(1,2,0,0) # the last two parameters are not to be estimated and fixed to be 1 as indicated in par0
#' fit1 <- fitCopulaOne(par0, patternpar=patternpar, dat=dat, lower=lower, upper=upper, copula_family="PPPP")
#' patternpar <- c(1,2,3,3) # the last two parameters are the same and are to be estimated
#' fit2 <- fitCopulaOne(par0, patternpar=patternpar, dat=dat, lower=lower, upper=upper, copula_family="PPPP")
#' }
fitCopulaOne <- function(par0, patternpar=seq(1,length(par0)), dat, flag=1, integration=F, opt="L-BFGS-B", se=F, lower=rep(0.1, length(unique(patternpar))), upper=rep(5, length(unique(patternpar))), trace=0, factr=1e9, printlevel=0, copula_family="PPPP")
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
      par0[which(patternpar!=0)] <- exp(par[patternpar[which(patternpar!=0)]])+0.01
      nllk_each_i <- try(parallel::parLapply(cl, seq_len(nco), seqRun, dat=dat, nco=nco, para=par0, flag=flag, integration=integration, copula_family=copula_family), silent = T)
      if(is(nllk_each_i,"try-error") | is(nllk_each_i,"warning")){return(20)}else{
        nllk_each_i <- unlist(nllk_each_i)
        out <- sum(nllk_each_i[is.finite(nllk_each_i)])
        return(out)
      }
    }
    if(se==T){hes <- T}else{hes <- F}
    fit <- nlm(obj, p=log(par0[match(setdiff(unique(patternpar),0), patternpar)]-0.01), print.level=printlevel, hessian = hes)
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
      par0[which(patternpar!=0)] <- par[patternpar[which(patternpar!=0)]]
      nllk_each_i <- try(parallel::parLapply(cl, seq_len(nco), seqRun, dat=dat, nco=nco, para=par0, flag=flag, integration=integration, copula_family=copula_family), silent = T)
      if(is(nllk_each_i,"try-error") | is(nllk_each_i,"warning")){return(20)}else{
        nllk_each_i <- unlist(nllk_each_i)
        out <- sum(nllk_each_i[is.finite(nllk_each_i)])
        return(out)
      }
    }
    if(se==T){hes <- T}else{hes <- F}
    fit <- optim(par=par0[match(setdiff(unique(patternpar),0), patternpar)], obj, method="L-BFGS-B", control=list(trace=trace, factr=factr), hessian=hes, lower=lower, upper=upper)
    parallel::stopCluster(cl)
    fit$par <- fit$par
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
    fit$warning <- "The second variable was transformed to 1 - v"
    fit$neg <- 1
  }else
  {
    fit$neg <- 0
  }
  fit
}

#' Contour plots of CopulaOne
#'
#' Contour plots based on either normal scores or uniform scores
#' @param para: vector of parameters (GGEE: a, b; PPPP: al, be, a, b)
#' @param marg: indicate normal scores or uniform scores (default: marg = "normal")
#' @param resolution: number of evaluations for each dimension (default is 30)
#' @param flag: used in the appell function (default is 1)
#' @param main: title of the plot, default: details of the plot
#' @keywords contour plot
#' @export
#' @examples
#' \dontrun{
#' plotCopulaOne(c(0.5, 1.8), copula_family="GGEE")
#' plotCopulaOne(c(0.5, 1.8), marg="uniform", resolution=20, copula_family="GGEE")
#' plotCopulaOne(c(0.5, 2.1,1,1), resolution=100, copula_family="PPPP")
#' }
plotCopulaOne <- function(para, marg="normal", drawlabels=T, flag=1, integration=F, resolution=30, copula_family="PPPP", main=NULL)
{
  zvec <- seq(-2.5, 2.5, length=resolution)
  f <- dnorm(zvec)
  Fvec <- pnorm(zvec)
  nn <- length(zvec)

  nco <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(nco)
  parallel::clusterEvalQ(cl, library(CopulaOne))

  Fmat <- cbind(rep(Fvec, each = nn), rep(Fvec, times = nn))

  if(copula_family=="GGEE" | copula_family=="PPPP")
  {
    denvec <- try(parallel::parLapply(cl, seq_len(nco), seqRun, dat=Fmat, nco=nco, para=para, flag=flag, integration=integration, copula_family=copula_family), silent = T)
  }else
  {
    cat("Please specify copula families from (GGEE, PPPP)", "\n")
    return(1)
  }
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
    if(length(main)==0)
    {
      main <- paste(copula_family, " (para=", paste(format(para), collapse=" "), ")", sep="")
    }
    contour(zvec,zvec, denmat, drawlabels=drawlabels, levels = pretty(rangeforplot, 10), labcex = 0.6, main=main, xlim = c(-3,3), ylim = c(-3,3))
  } else if(marg == "uniform")
  {
    denmat <- matrix(denvec, nn, nn)
    rangeforplot <- denvec
    if(length(main)==0)
    {
      main <- paste(copula_family, " (para=", paste(format(para), collapse=" "), ")", sep="")
    }
    contour(Fvec, Fvec, denmat, drawlabels=drawlabels, levels = pretty(rangeforplot, 20), labcex = 0.6, main=main, xlim = c(0,1), ylim = c(0,1))
  }
  invisible(denvec)
}
