# copulaone.R

.copula_worker_count <- function(n, cores = parallel::detectCores()) {
  if (length(n) != 1L || !is.finite(n) || n < 1 || n != floor(n))
    stop("At least one observation is required.")
  if (length(cores) != 1L || !is.finite(cores)) cores <- 1L
  as.integer(min(n, max(1, cores - 1)))
}

# A finite penalty is required by L-BFGS-B; never omit failed observations.
.complete_copula_nllk <- function(parts, n) {
  if (inherits(parts, "try-error") || inherits(parts, "warning")) return(1e100)
  values <- unlist(parts, use.names = FALSE)
  if (!is.numeric(values) || length(values) != n || any(!is.finite(values)))
    return(1e100)
  value <- sum(values)
  if (is.finite(value)) value else 1e100
}

# used for parallel computing
seqRun <- function(i, dat, nco, para, flag=1, integration=F, copula_family="PPPP")
{
  # dataset is divided in nco parts, and the ith part is working
  dat <- as.matrix(dat)
  N <- dim(dat)[1]
  if (length(nco) != 1L || !is.finite(nco) || nco < 1 || nco != floor(nco) ||
      length(i) != 1L || !is.finite(i) || i < 1 || i > nco || i != floor(i))
    stop("i and nco must identify a valid worker.")
  first <- floor((i - 1) * N / nco) + 1L
  last <- floor(i * N / nco)
  if (first > last) return(numeric(0))
  dat_i <- dat[seq.int(first, last), , drop = FALSE]
  if(copula_family=="PPPP")
  {
    den <- dPPPP_COP(u=dat_i[,1], v=dat_i[,2], al=para[1], be=para[2], a=para[3], b=para[4])
  }else if(copula_family=="GGEE")
  {
    den <- dGGEE_COP(u=dat_i[,1], v=dat_i[,2], al=para[1], be=para[2], flag=flag, integration=integration)
  }else
  {
    stop("Please specify copula families from (GGEE, PPPP).")
  }
  if (!is.numeric(den) || length(den) != nrow(dat_i) ||
      any(!is.finite(den) | den <= 0))
    stop("Every observation must have a finite positive copula density.")
  return(-log(den))
}

# FRA1 fitting uses the vectorized log-density from fullrange2.R.
# No worker cluster is needed.
.fitFRA1_CopulaOne <- function(par0, patternpar, dat, opt, se,
                                lower, upper, trace, factr, printlevel) {
  if (!opt %in% c("L-BFGS-B", "nlm"))
    stop("For FRA1, opt must be 'L-BFGS-B' or 'nlm'.")
  dat <- as.matrix(dat)
  if (!is.numeric(dat) || ncol(dat) != 2L || nrow(dat) < 2L ||
      any(!is.finite(dat)) || any(dat <= 0 | dat >= 1))
    stop("dat must have two numeric columns of finite scores strictly in (0,1).")
  if (!is.numeric(par0) || length(par0) != 2L ||
      any(!is.finite(par0)) || any(par0 < -1 | par0 >= 1))
    stop("FRA1 par0 must be c(eta, theta), with each value in [-1,1).")
  if (!is.numeric(patternpar) || length(patternpar) != 2L ||
      any(!is.finite(patternpar)) || any(patternpar < 0 | patternpar != floor(patternpar)))
    stop("patternpar must contain two nonnegative integer group labels.")
  groups <- sort(unique(patternpar[patternpar > 0]))
  k <- length(groups)
  index <- match(patternpar, groups)
  free <- patternpar > 0
  par_names <- c("eta", "theta")
  free_names <- vapply(groups, function(g)
    paste(par_names[patternpar == g], collapse = "="), character(1))
  start <- par0[match(groups, patternpar)]
  bound <- function(x, name) {
    if (k == 0L) return(numeric(0))
    if (length(x) == 1L) x <- rep(x, k)
    if (length(x) != k || any(!is.finite(x)))
      stop(name, " must have one value or one value per free parameter group.")
    as.numeric(x)
  }
  lower <- bound(lower, "lower")
  upper <- bound(upper, "upper")
  if (any(lower < -1 | upper >= 1 | lower > upper))
    stop("FRA1 bounds must satisfy -1 <= lower <= upper < 1.")
  if (any(start < lower | start > upper))
    stop("Initial free parameters must lie within lower and upper.")
  unpack <- function(p) {
    full <- par0
    full[free] <- p[index[free]]
    setNames(full, par_names)
  }
  ken <- suppressWarnings(cor(dat[,1], dat[,2], method = "kendall"))
  if (!is.finite(ken)) stop("Kendall correlation is undefined; check constant margins.")
  neg <- as.integer(ken < 0)
  if (neg) dat[,2] <- 1 - dat[,2]
  objective <- function(p) {
    full <- unpack(p)
    ll <- tryCatch(logdFRA1_COP(dat[,1], dat[,2],
                                 eta = full[1], theta = full[2]), error = function(e) NA_real_)
    if (length(ll) != nrow(dat) || any(!is.finite(ll))) return(1e100)
    value <- -sum(ll)
    if (is.finite(value)) value else 1e100
  }
  # Separate negative, zero and positive intervals for each free group.
  # This includes the transition axes and supports tied eta/theta parameters.
  parts <- lapply(seq_len(k), function(j) {
    lo <- lower[j]; hi <- upper[j]; ans <- list()
    if (lo < 0) ans[[length(ans) + 1L]] <- c(lo, min(hi, 0))
    if (lo <= 0 && hi >= 0) ans[[length(ans) + 1L]] <- c(0, 0)
    if (hi > 0) ans[[length(ans) + 1L]] <- c(max(lo, 0), hi)
    ans
  })
  grid <- if (k) expand.grid(lapply(parts, seq_along)) else data.frame(dummy = 1L)
  candidates <- list()
  for (row in seq_len(nrow(grid))) {
    lo <- hi <- numeric(k)
    for (j in seq_len(k)) {
      interval <- parts[[j]][[grid[row,j]]]
      lo[j] <- interval[1]; hi[j] <- interval[2]
    }
    active <- which(lo < hi)
    assemble <- function(x) { p <- lo; p[active] <- x; p }
    if (!length(active)) {
      candidates[[length(candidates) + 1L]] <- list(
        par = lo, value = objective(lo), convergence = 0L,
        message = "Fixed parameters/transition candidate")
      next
    }
    branch_obj <- function(x) objective(assemble(x))
    starts <- list(pmin(pmax(start[active], lo[active]), hi[active]),
                   (lo[active] + hi[active]) / 2)
    for (s in starts) {
      result <- tryCatch({
        if (opt == "L-BFGS-B") {
          z <- optim(s, branch_obj, method = "L-BFGS-B",
                     lower = lo[active], upper = hi[active],
                     control = list(trace = trace, factr = factr))
          list(par = assemble(z$par), value = z$value,
               convergence = z$convergence, message = z$message)
        } else {
          # A logistic transform enforces each branch's finite bounds.
          width <- hi[active] - lo[active]
          decode <- function(z) lo[active] + width * plogis(z)
          z0 <- qlogis(pmin(pmax((s - lo[active]) / width, 1e-6), 1 - 1e-6))
          z <- nlm(function(x) branch_obj(decode(x)), p = z0,
                   print.level = printlevel, hessian = FALSE)
          list(par = assemble(decode(z$estimate)), value = z$minimum,
               convergence = if (z$code <= 2L) 0L else z$code,
               message = paste("nlm code", z$code), code = z$code,
               iterations = z$iterations)
        }
      }, error = function(e) list(par = assemble(s), value = Inf,
                                  convergence = 99L, message = conditionMessage(e)))
      candidates[[length(candidates) + 1L]] <- result
    }
  }
  values <- vapply(candidates, function(x) x$value, numeric(1))
  best <- which.min(values)
  if (!length(best) || !is.finite(values[best]) || values[best] >= 1e100)
    stop("FRA1 fitting failed: no candidate has a finite complete log-likelihood.")
  fit <- candidates[[best]]
  names(fit$par) <- free_names
  fit$fullpar <- unpack(fit$par)
  fit$loglik <- -fit$value
  fit$copula_family <- "FRA1"
  fit$patternpar <- patternpar
  fit$neg <- neg
  fit$branch_results <- data.frame(
    nllk = values,
    convergence = vapply(candidates, function(x) x$convergence, integer(1)))
  if (opt == "nlm") {
    fit$estimate <- fit$par
    fit$minimum <- fit$value
    if (is.null(fit$code)) fit$code <- 1L
  }
  if (fit$convergence != 0L)
    warning("The best FRA1 candidate did not report convergence; inspect branch_results.")
  if (neg) {
    fit$warning <- "The second variable was transformed to 1 - v"
    cat(fit$warning, "\n")
  }
  if (se) {
    fit$err <- setNames(rep(NA_real_, k), free_names)
    # Central Hessians are inappropriate at transitions or active bounds.
    near_edge <- any(abs(fit$par) <= 0.005 |
                       fit$par - lower <= 0.005 | upper - fit$par <= 0.005)
    if (k && !near_edge && fit$convergence == 0L) {
      h <- tryCatch(optimHess(fit$par, objective), error = function(e) NULL)
      if (!is.null(h) && all(is.finite(h)) &&
          min(eigen(h, symmetric = TRUE, only.values = TRUE)$values) > 0) {
        fit$hessian <- h
        fit$err <- setNames(sqrt(diag(solve(h))), free_names)
      } else warning("FRA1 Hessian is unavailable or not positive definite; SEs are NA.")
    } else if (k) {
      warning("FRA1 Hessian SEs are NA near transitions/bounds or after nonconvergence; use profile likelihood or bootstrap.")
    }
  }
  fit
}
#' Model fitting with CopulaOne
#'
#' Fit bivariate data with the full-range tail dependence copula based on maximum likelihood
#' @param par0 Initial parameters: eta, theta for FRA1; al, be for GGEE; al, be, a, b for PPPP.
#' @param patternpar Nonnegative integer group labels, one per parameter. Zero fixes a parameter at par0; equal positive labels tie parameters. Groups are ordered by label.
#' @param dat Two numeric columns of finite scores strictly in (0,1), with at least two observations.
#' @param flag Method flag passed to appell::appellf1().
#' @param integration Logical; use numerical integration instead of Appell's function.
#' @param opt Optimization method: L-BFGS-B or nlm. Both enforce the supplied bounds.
#' @param se Logical; request Hessian standard errors. Unavailable estimates are NA with a warning.
#' @param lower Lower bound, one scalar or one value per free group. Positive for GGEE/PPPP; at least -1 for FRA1.
#' @param upper Upper bound, one scalar or one value per free group. Must be less than 1 for FRA1.
#' @param trace Trace level for optim().
#' @param factr Convergence tolerance multiplier for L-BFGS-B.
#' @param printlevel Print level for nlm().
#' @param copula_family Copula family: FRA1 (default), GGEE, or PPPP.
#' @param workers Positive worker count for GGEE/PPPP; defaults to option CopulaOne.workers or 1. FRA1 is evaluated serially.
#' @returns A list with free parameters par, fullpar, value (negative log-likelihood), loglik, convergence, neg, and optional err. For nlm, estimate and minimum are also returned. FRA1 additionally returns branch_results.
#' @export
#' @keywords distribution
#' @examples
#' \dontrun{
#' uv <- rFRA1_COP(500, eta=0.3, theta=0.4, seed=123)
#' fit <- fitCopulaOne(dat=uv) # FRA1, with default initial values
#' fit$fullpar # eta, theta
#' fit <- fitCopulaOne(c(0.3, 0.4), dat=uv, copula_family="FRA1")
#' data("euro0306")
#' dat <- uscore(euro0306[,c(2,3)])[1:100,]
#' par0 <- c(0.3, 0.3)
#' fit <- fitCopulaOne(par0, dat=dat, copula_family="GGEE")
#' par0 <- c(0.3, 0.3, 1, 1)
#' lower <- 0.1
#' upper <- 5
#' par0 <- c(0.3, 0.3, 1, 1)
#' patternpar <- c(1,2,0,0) # the last two parameters are not to be estimated and fixed to be 1 as indicated in par0
#' fit1 <- fitCopulaOne(par0, patternpar=patternpar, dat=dat, lower=lower, upper=upper, copula_family="PPPP")
#' patternpar <- c(1,2,3,3) # the last two parameters are the same and are to be estimated
#' fit2 <- fitCopulaOne(par0, patternpar=patternpar, dat=dat, lower=lower, upper=upper, copula_family="PPPP")
#' }
fitCopulaOne <- function(par0=if(copula_family=="PPPP") c(.2,.2,1,1) else c(.2,.2),
  patternpar=seq_along(par0), dat, flag=1, integration=FALSE, opt="L-BFGS-B", se=FALSE,
  lower=if(copula_family=="FRA1") -1 else .1,
  upper=if(copula_family=="FRA1") .95 else 5,
  trace=0, factr=1e9, printlevel=0, copula_family="FRA1",
  workers=getOption("CopulaOne.workers",1L)) {
  copula_family<-match.arg(copula_family,c("FRA1","GGEE","PPPP"))
  opt<-match.arg(opt,c("L-BFGS-B","nlm"))
  dat<-as.matrix(dat)
  if(!is.numeric(dat) || ncol(dat)!=2L || nrow(dat)<2L ||
      any(!is.finite(dat) | dat<=0 | dat>=1))
    stop("dat must have two numeric columns of finite scores strictly in (0,1).")
  if(length(workers)!=1L || !is.numeric(workers) || !is.finite(workers) ||
      workers<1 || workers!=floor(workers)) stop("workers must be a positive integer.")
  if(!is.logical(se) || length(se)!=1L || is.na(se) ||
      !is.logical(integration) || length(integration)!=1L || is.na(integration))
    stop("se and integration must be TRUE or FALSE.")
  if(copula_family=="FRA1")
    return(.fitFRA1_CopulaOne(par0,patternpar,dat,opt,se,lower,upper,trace,factr,printlevel))
  .fit_positive_copula(par0,patternpar,dat,opt,se,lower,upper,trace,factr,printlevel,
                       copula_family,flag,integration,as.integer(workers))
}
#' Contour plots of CopulaOne
#'
#' Contour plots based on either normal scores or uniform scores
#' @param para Parameters: eta, theta for FRA1; al, be for GGEE; al, be, a, b for PPPP.
#' @param marg Plot margins: normal or uniform.
#' @param drawlabels Logical; label contour lines.
#' @param flag Method flag passed to appell::appellf1().
#' @param integration Logical; use numerical integration instead of Appell's function.
#' @param resolution Number of grid evaluations per dimension, an integer of at least 2.
#' @param copula_family Copula family: FRA1 (default), GGEE, or PPPP.
#' @param main Optional plot title.
#' @returns The evaluated density vector, invisibly; draws a contour plot.
#' @export
#' @keywords distribution
#' @examples
#' \dontrun{
#' plotCopulaOne(c(0.6, -0.4)) # FRA1: eta=0.6, theta=-0.4
#' plotCopulaOne(c(0.6, -0.4), marg="uniform")
#' plotCopulaOne(c(0.5, 1.8), copula_family="GGEE")
#' plotCopulaOne(c(0.5, 1.8), marg="uniform", resolution=20, copula_family="GGEE")
#' plotCopulaOne(c(0.5, 2.1,1,1), resolution=100, copula_family="PPPP")
#' }
plotCopulaOne <- function(para,marg="normal",drawlabels=TRUE,flag=1,integration=FALSE,
                           resolution=30,copula_family="FRA1",main=NULL) {
  copula_family<-match.arg(copula_family,c("FRA1","GGEE","PPPP"))
  marg<-match.arg(marg,c("normal","uniform"))
  if(!is.numeric(resolution) || length(resolution)!=1L || !is.finite(resolution) ||
      resolution<2 || resolution!=floor(resolution)) stop("resolution must be an integer of at least 2.")
  count<-if(copula_family=="PPPP") 4L else 2L
  if(!is.numeric(para) || length(para)!=count || any(!is.finite(para)))
    stop("para has an invalid length or contains non-finite values.")
  if(copula_family=="FRA1") {
    check_parameter_FRA1(para[1],"eta"); check_parameter_FRA1(para[2],"theta")
  } else if(any(para<=0)) stop("GGEE/PPPP parameters must be positive.")
  zvec<-seq(-2.5,2.5,length.out=resolution)
  Fvec<-pnorm(zvec); nn<-length(zvec)
  Fmat<-cbind(rep(Fvec,each=nn),rep(Fvec,times=nn))
  if(copula_family=="FRA1") denvec<-dFRA1_COP(Fmat[,1],Fmat[,2],para[1],para[2]) else {
    workers<-getOption("CopulaOne.workers",1L)
    if(length(workers)!=1 || !is.finite(workers) || workers<1 || workers!=floor(workers))
      stop("CopulaOne.workers must be a positive integer.")
    nco<-min(workers,nrow(Fmat))
    if(nco==1L) denvec<-exp(-seqRun(1,Fmat,1,para,flag,integration,copula_family)) else {
      cl<-parallel::makeCluster(nco)
      on.exit(parallel::stopCluster(cl),add=TRUE)
      parallel::clusterCall(cl,function(paths) {.libPaths(paths); loadNamespace("CopulaOne");NULL},.libPaths())
      denvec<-exp(-unlist(parallel::parLapply(cl,seq_len(nco),seqRun,dat=Fmat,nco=nco,
                         para=para,flag=flag,integration=integration,copula_family=copula_family)))
    }
  }
  if(length(denvec)!=nrow(Fmat) || any(!is.finite(denvec))) stop("Density failed on the plotting grid.")
  if(is.null(main)) main<-paste0(copula_family," (",paste(format(para),collapse=", "),")")
  if(marg=="normal") denvec<-denvec*rep(dnorm(zvec),each=nn)*rep(dnorm(zvec),times=nn)
  axis<-if(marg=="normal") zvec else Fvec
  contour(axis,axis,matrix(denvec,nn,nn),drawlabels=drawlabels,main=main)
  invisible(denvec)
}
