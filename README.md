# CopulaOne - an R package for full-range tail dependence copulas


The R package *CopulaOne* implements functions for bivariate full-range tail dependence copulas.

Bivariate copulas have been widely used either in modeling bivariate dependence structures or building multivariate dependence models such as Vine copulas and factor copulas. In the literature, there are many parametric bivariate copula families, but they often have specific dependence patterns, which limit their use in real applications. The R package *CopulaOne* aims at implementing a collection of very flexible bivariate copulas that are parsimonious and very flexible. The copulas implemented in *CopulaOne* should be able to account for most bivariate dependence patterns by a single copula, and this is also why we name the package as *CopulaOne*. Compared to those existing bivariate parametric copula families, the main merit of the bivariate copulas implemented here is that, they can account for full-range tail dependence in both upper and lower tails. The package is under active development, and the following copulas have been implemented: GGEE, PPPP, FRA1. The following figures show how flexible these copulas are. Please refer to [Hua (2017)](https://doi.org/10.1016/j.insmatheco.2017.01.003) for more details about the GGEE copula, [Su and Hua (2017)](https://doi.org/10.1016/j.insmatheco.2017.08.009) for the details about the PPPP copula, and [Hua (2026)](https://arxiv.org/abs/2609.18742) for the details about the FRA1 copula.

|                       |                       |                       |
| --------------------- | --------------------- | --------------------- |
| <img src="https://github.com/larryleihua/CopulaOne/blob/master/inst/extdata/FRA1.gif" width="200" height="224" />  | <img src="https://github.com/larryleihua/CopulaOne/blob/master/inst/extdata/GGEE.gif" width="200" height="224" />  |  <img src="https://github.com/larryleihua/CopulaOne/blob/master/inst/extdata/PPPP.gif" width="200" height="224" /> |


CopulaOne implements bivariate copulas with flexible upper and lower tail dependence.
The default family is **FRA1**; **GGEE** and **PPPP** are also available.

## Installation

Install `hypergeo`, `appell`, and `cubature` in the library used by your R session,
then install this source package:

```r
install.packages(c("hypergeo", "appell", "cubature"))
install.packages("remotes")
remotes::install_github("larryleihua/CopulaOne", force=T)

devtools::install_github("larryleihua/CopulaOne", force=T) # alternative method
```

If `appell` is unavailable from your repository, its source archive is an alternative:

```r
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/appell/appell_0.0-4.tar.gz",
  repos = NULL, type = "source"
)
```

- If there are some issues from the above codes, you may need to install the following R packages first: _appell_ and _hypergeo_.

_appell_ can be installed by the following R codes:
```{r, eval=FALSE}
install.packages("appell_0.0-4.tar.gz", repos = NULL, type = "source")
```
where appell_0.0-4.tar.gz can be downloaded from the following website: (accessible Jan. 31, 2025)
https://cran.r-project.org/src/contrib/Archive/appell/

1. If you use Windows OS, then you will need to install Rtools in advance.
   
2. If you use MacOS, then try the following steps on terminal to install gfortran and its paths:
     1. Install Homebrew
        ```
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
        ```
     1. Install gfortran (Fortran Compiler) via Homebrew
        ```
        brew install gcc
        ```
     1. Create ~/.R/Makevars, and add the following to the file (find the paths on your computer first):
        ```
        FC = /opt/homebrew/bin/gfortran
        F77 = /opt/homebrew/bin/gfortran
        LDFLAGS += -L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/14
        FLIBS = -lgfortran
        ```
   
        - use following to find your path '/opt/homebrew/bin/gfortran'
          ```
          which gfortran
          ```          
        - use following to find path '/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/14'
          ```
          find /opt/homebrew -name "libgfortran*"
          ```
     1. On R, use following to install _appell_
        ```
        install.packages("Pathto/appell_0.0-4.tar.gz", repos = NULL, type = "source")
        ```

_hypergeo_ can be installed easily:
```{r, eval=FALSE}
install.packages("hypergeo", dependencies = T)
```

Installing compiled dependencies from source requires the compiler toolchain for
that R installation (Rtools on Windows). Check `R.version.string` and `.libPaths()`
when an installed dependency cannot be found. The package requires R >= 4.3;
its regression checks have been exercised with R 4.5.3.

## FRA1

Parameters are ordered `eta` (lower tail), then `theta` (upper tail), each in
`[-1,1)`. Independence is `eta = theta = -1`.

```r
library(CopulaOne)
uv <- rFRA1_COP(800, eta = 0.3, theta = 0.4, seed = 123)
pFRA1_COP(0.3, 0.4, eta = 0.3, theta = 0.4)
dFRA1_COP(0.3, 0.4, eta = 0.3, theta = 0.4)
logdFRA1_COP(0.3, 0.4, eta = 0.3, theta = 0.4)
C2FRA1_COP(0.3, 0.4, eta = 0.3, theta = 0.4)
C2invFRA1_COP(0.5, 0.4, eta = 0.3, theta = 0.4)
tailFRA1_COP(eta = 0.3, theta = 0.4)
dependenceFRA1_COP(eta = 0.3, theta = 0.4)
fit <- fitCopulaOne(dat = uv)
fit$fullpar
plotCopulaOne(fit$fullpar)
```

## GGEE and PPPP

GGEE uses positive parameters `al, be`; PPPP uses positive parameters
`al, be, a, b`. The distribution and density functions use the same parameter order.

```r
rGGEE_COP(10, al = 0.5, be = 0.8, seed = 1)
rPPPP_COP(10, al = 0.5, be = 0.8, a = 1, b = 1, seed = 1)
dGGEE_COP(0.2, 0.3, al = 1.2, be = 0.5)
pPPPP_COP(0.2, 0.3, al = 1.2, be = 0.5, a = 1, b = 1)
```

## Fitting and uniform scores

```r
data(euro0306)
x <- euro0306[, c(2, 3)]
x <- x[complete.cases(x), ]
uv <- uscore(x)
fit_PPPP <- fitCopulaOne(c(0.3, 0.3, 1, 1), dat = uv,
                         patternpar = c(1, 2, 0, 0), copula_family = "PPPP")
fit_FRA1 <- fitCopulaOne(dat=uv)
```

GGEE/PPPP fitting defaults to one worker. Use `workers = 2` or set
`options(CopulaOne.workers = 2)` to enable parallel evaluation. FRA1 is vectorized
and runs serially. Worker processes inherit the current library paths.

## Development checks

Run from the repository root with the intended R installation and library:

```sh
Rscript scripts/check.R
```

The script builds a source archive and runs `R CMD check --no-manual` in a
temporary directory, including the regression tests. A Windows CI workflow
runs the same checks. Generate help files after changing roxygen comments with
`roxygen2::roxygenise()`.

## References

GGEE: [Hua (2017)](https://doi.org/10.1016/j.insmatheco.2017.01.003).
PPPP: [Su and Hua (2017)](https://doi.org/10.1016/j.insmatheco.2017.08.009).
FRA1: [Hua (2026)](https://arxiv.org/abs/2609.18742).

## Citation, please use the following bibtex for citation

```
@misc{Hua2026,
  author = {Lei Hua},
  title  = {Copula{O}ne - an {R} package for full-range tail dependence copulas},
  year   = {2026},
  howpublished = "\url{https://github.com/larryleihua/CopulaOne}"
}
```
