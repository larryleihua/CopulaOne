.fit_positive_copula <- function(par0, patternpar, dat, opt, se, lower, upper,
                                  trace, factr, printlevel, copula_family,
                                  flag, integration, workers) {
  count <- if (copula_family=="GGEE") 2L else 4L
  if (!is.numeric(par0) || length(par0)!=count || any(!is.finite(par0) | par0<=0))
    stop("par0 must contain ",count," finite positive parameters.")
  if (!is.numeric(patternpar) || length(patternpar)!=count ||
      any(!is.finite(patternpar) | patternpar<0 | patternpar!=floor(patternpar)))
    stop("patternpar must have one nonnegative integer label per parameter.")
  groups <- sort(unique(patternpar[patternpar>0]))
  k <- length(groups); free <- patternpar>0
  index <- match(patternpar,groups)
  start <- par0[match(groups,patternpar)]
  par_names <- if(count==2L) c("al","be") else c("al","be","a","b")
  free_names <- vapply(groups,function(g) paste(par_names[patternpar==g],collapse="="),character(1))
  bound <- function(x,name) {
    if (!k) return(numeric(0))
    if(length(x)==1L) x<-rep(x,k)
    if(!is.numeric(x) || length(x)!=k || any(!is.finite(x) | x<=0))
      stop(name," must contain a positive finite value per free group, or one scalar.")
    x
  }
  lower<-bound(lower,"lower"); upper<-bound(upper,"upper")
  if(any(lower>upper | start<lower | start>upper)) stop("Invalid bounds or initial values outside bounds.")
  unpack<-function(p) {
    full<-par0; full[free]<-p[index[free]]; setNames(full,par_names)
  }
  ken<-cor(dat[,1],dat[,2],method="kendall")
  if(!is.finite(ken)) stop("Kendall correlation is undefined; check constant margins.")
  neg<-as.integer(ken<0)
  if(neg) dat[,2]<-1-dat[,2]
  nco<-min(workers,nrow(dat))
  cl<-NULL
  if(nco>1L) {
    cl<-parallel::makeCluster(nco)
    on.exit(parallel::stopCluster(cl),add=TRUE)
    parallel::clusterCall(cl,function(paths) {
      .libPaths(paths); loadNamespace("CopulaOne"); NULL
    },.libPaths())
  }
  objective<-function(p) {
    parts<-tryCatch({
      if(is.null(cl)) list(seqRun(1L,dat,1L,unpack(p),flag,integration,copula_family)) else
        parallel::parLapply(cl,seq_len(nco),seqRun,dat=dat,nco=nco,
                            para=unpack(p),flag=flag,integration=integration,
                            copula_family=copula_family)
    },error=function(e) structure(conditionMessage(e),class="try-error"))
    .complete_copula_nllk(parts,nrow(dat))
  }
  active<-which(lower<upper)
  assemble<-function(x) {p<-lower; p[active]<-x; p}
  obj<-function(x) objective(assemble(x))
  if(!length(active)) {
    fit<-list(par=lower,value=objective(lower),convergence=0L)
  } else if(opt=="L-BFGS-B") {
    fit<-optim(start[active],obj,method=opt,lower=lower[active],upper=upper[active],
               control=list(trace=trace,factr=factr))
    fit$par<-assemble(fit$par)
  } else {
    width<-upper[active]-lower[active]
    decode<-function(z) lower[active]+width*plogis(z)
    z0<-qlogis(pmin(pmax((start[active]-lower[active])/width,1e-6),1-1e-6))
    fit<-nlm(function(z) obj(decode(z)),p=z0,print.level=printlevel)
    fit$par<-assemble(decode(fit$estimate)); fit$value<-fit$minimum
    fit$convergence<-if(fit$code<=2L) 0L else fit$code
  }
  if(!is.finite(fit$value) || fit$value>=1e100)
    stop("Fitting failed: no finite complete log-likelihood at the returned parameters.")
  names(fit$par)<-free_names
  fit$fullpar<-unpack(fit$par); fit$loglik<--fit$value
  fit$copula_family<-copula_family; fit$patternpar<-patternpar; fit$neg<-neg
  if(opt=="nlm") {fit$estimate<-fit$par; fit$minimum<-fit$value}
  if(fit$convergence!=0L) warning("Optimizer did not converge; inspect the fit result.")
  if(neg) fit$warning<-"The second variable was transformed to 1 - v"
  if(se) {
    fit$err<-setNames(rep(NA_real_,k),free_names)
    interior<-k>0L && all(fit$par-lower>.005 & upper-fit$par>.005)
    if(interior && fit$convergence==0L) {
      h<-tryCatch(optimHess(fit$par,objective),error=function(e) NULL)
      inv<-tryCatch({
        if(is.null(h) || any(!is.finite(h)) || min(eigen(h,symmetric=TRUE,only.values=TRUE)$values)<=0)
          stop("Invalid Hessian")
        solve(h)
      },error=function(e) NULL)
      if(!is.null(inv)) {fit$hessian<-h; fit$err<-setNames(sqrt(diag(inv)),free_names)}
    }
    if(k && anyNA(fit$err)) warning("Hessian SEs are unavailable at bounds, after nonconvergence, or for a singular Hessian.")
  }
  fit
}
