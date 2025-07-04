% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fullrange2.R
\name{pGGEE_COP}
\alias{pGGEE_COP}
\title{Joint CDF of the GGEE copula model}
\usage{
pGGEE_COP(u, v, al, be, flag = 1, integration = F, maxit = 1e+05)
}
\arguments{
\item{u, v:}{data input.}

\item{al, be:}{parameters.}

\item{flag:}{flag used in the appell::appellf1() function.}
}
\description{
Joint CDF of the GGEE copula model
}
\examples{
pGGEE_COP(0.9, 0.3, 1, 1)
pGGEE_COP(0.9, 0.3, 1, 1, integration=TRUE)
}
\keyword{CDF}
\keyword{Joint}
