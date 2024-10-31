# CopulaOne - an R package for full-range tail dependence copulas

Warning: This is the development version of the package, please use with caution!

The R package *CopulaOne* implements functions for bivariate copulas that must satisfy the following two properties:
* It can account for full-range tail dependence for both upper and lower tails.
* It can account for both reflection symmetry and asymmetry between upper and lower tails.

Bivariate copulas have been widely used either in modeling bivariate dependence structures or building multivariate dependence models such as Vine copulas and factor copulas. In the literature, there are numerous parametric bivariate copula families. It is often very time consuming to select copula families from many different candidate copula families. The R package *CopulaOne* aims at implementing a collection of very flexible bivariate copulas that are parsimonious and very flexible. The copulas implemented in *CopulaOne* should be able to account for most bivariate dependence patterns by a single copula, and this is also why we name the package as *CopulaOne*. Compared to those existing bivariate parametric copula families, the main merit of the bivariate copulas implemented here is that, they can account for full-range tail dependence in both upper and lower tails, and the upper and lower tails can be either reflection symmetric or asymmetric.
The package is under active development, and the following copulas have been implemented: GGEE, PPPP. The following figures show how flexible the GGEE and PPPP copulas are. Please refer to [Hua (2017)](https://doi.org/10.1016/j.insmatheco.2017.01.003) for more details about the GGEE copula, and [Su and Hua (2017)](https://doi.org/10.1016/j.insmatheco.2017.08.009) for the details about the PPPP copula.

|                           |                           | 
| ------------------------- | ------------------------- |
| <img src="https://github.com/larryleihua/CopulaOne/blob/master/data/GGEE.gif" width="250" height="280" />  |  <img src="https://github.com/larryleihua/CopulaOne/blob/master/data/PPPP.gif" width="250" height="280" /> |

## Installation
- The R package CopulaOne can be easily installed from github by the following two lines.
```{r, eval=FALSE}
library(devtools)
install_github("larryleihua/CopulaOne", force=T)
```
- If there are some issues from the above codes, you may need to install the following R packages first: _appell_ and _hypergeo_.

_appell_ can be installed by the following R codes:
```{r, eval=FALSE}
install.packages("appell_0.0-4.tar.gz", repos = NULL, type = "source")
```
where appell_0.0-4.tar.gz can be downloaded from the following website: (accessible Oct. 1, 2023)
https://cran.r-project.org/src/contrib/Archive/appell/
If you use Windows OS, then you will need to install Rtools in advance.


_hypergeo_ can be installed easily:
```{r, eval=FALSE}
install.packages("hypergeo", dependencies = T)
```

## Basic functions
- Naming rules: The name *GGEE_COP* is used for the two-parameter copula that are based on Gamma-Gamma-Exponential-Exponential mixtures. The name *CopulaOne* is used as a unified platform for implementing various functions that can be used as coherent as possible. For other copulas, replace GGEE by the corresponding names, such as PPPP.

- Simulation based on the copula can be done as follows:
```{r}
library(CopulaOne)
UU <- rGGEE_COP(10, a=0.5, b=0.8)
```

- Joint density, cdf functions can be derived as
```{r}
den <- dGGEE_COP(0.2, 0.3, 1.2, 0.5)
cdf <- pGGEE_COP(0.2, 0.3, 1.2, 0.5)
cat("The copula density and cdf are:", den, cdf, "\n")
```

- Contour plots can be plotted directly for a given copula as follows.
```{r fig.width=11, fig.height=6}
layout(matrix(c(1,2),1,2))
plotCopulaOne(c(1.2, 0.5), marg = "normal", copula_family = "GGEE")
plotCopulaOne(c(1.2, 0.5), marg = "uniform", copula_family = "GGEE")
```

- Kendall's tau and Spearman's rho of the GGEE copula can be evaluated by
```{r}
tauGGEE_COP(a=0.7, b=0.4)
sprGGEE_COP(a=0.7, b=0.4)
```

- Upper extreme value copula of the GGEE copula can also be evaluated with the following joint cdf and pdf functions. Lower extreme value copula should be the same if one exchanges a and b.
```{r}
pUEV_GGEE_COP(0.3, 0.4, b=1.2)
dUEV_GGEE_COP(0.3, 0.4, b=1.2)
```

## Model fitting
- An example of fitting dependence between exchange rates returns by the GGEE copula [Warning: this step can be slow on your computer!]
```{r}
data("euro0306")
dat <- uscore(euro0306[,c(2,3)])[1:50,]
par <- c(0.3, 0.3)
fit <- fitCopulaOne(par, dat=dat, copula_family = "GGEE")
```
- An example of fitting dependence between exchange rates returns by the PPPP copula 
```{r}
data("euro0306")
dat <- uscore(euro0306[,c(2,3)])[1:50,]
par0 <- c(0.3,0.3,1,1)
patternpar <- c(1,2,0,0)
fit1 <- fitCopulaOne(par0, patternpar=patternpar, dat=dat, se=F, copula_family = "PPPP")
```

## Known issues (i.e., to-do list)
- [solved] rGGEE_COP has issues when al and/or be are too small, say, 0.01, and rgamma() will generate vary small values so hypergeo::hypergeo will generate lots of boundary values 2.0, making the copula not working.


## Citation, please use the following bibtex for citation

```
@misc{Hua2018,
  author = {Lei Hua},
  title  = {Copula{O}ne - an {R} package for full-range tail dependence copulas},
  year   = {2018},
  howpublished = "\url{https://github.com/larryleihua/CopulaOne}"
}
```
