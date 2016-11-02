# CopulaOne - an R package for a very flexible bivariate parametric copula

The R package *CopulaOne* implements functions for a bivariate copula that is very flexible and parsimonious. 

Bivariate copulas have been widely used either in modeling bivariate dependence structures or building multivariate dependence models such as Vine copulas and factor copulas. In the literature, there are numerous parametric bivariate copula families. It is often very time consuming to select copula families from many different candidate copula families. The R package CopulaOne aims at implementing a very flexible bivariate copula that is parsimonious but very flexible and has the capacity of capturing most bivariate dependence patterns. Compared to those existing bivariate parametric copula families, the main merits of the bivariate copula is that, it can account for full-range tail dependence in both upper and lower tails, and the upper and lower tails can be either symmetric or asymmetric. 

The following shows how flexible the GGEE copula based on which the R package *CopulaOne* was built on. Refer to [Hua (2016)](http://niu.edu/leihua/papers/Manuscript-Full-Range-HUA.pdf) for more details about the GGEE copula.

![Animation of the GGEE copula](http://niu.edu/leihua/_images/copula.gif)

<img src="http://niu.edu/leihua/_images/copula.gif" width="200" height="200" />

## Installation
- The R package CopulaOne can be easily installed from github by the following two lines.
```{r, eval=FALSE}
library(devtools)
install_github("larryleihua/CopulaOne")
```

## Basic functions
- Naming rules: The name *GGEE_COP* is used for a specific two-parameter copula that are of full-range tail dependence in both upper and lower tails. The name *CopulaOne* is used as a unified platform for implementing various functions that can be used as coherent as possible.

- Simulation based on the copula can be done as follows:
```{r}
library(CopulaOne)
UU = rGGEE_COP(10, a=0.5, b=0.8)
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
plotCopulaOne(1.2, 0.5, marg = "normal")
plotCopulaOne(1.2, 0.5, marg = "uniform")
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
- An example of fitting dependence between exchange rates returns
```{r}
data("euro0306")
dat <- uscore(euro0306[,c(2,3)])[1:50,]
par <- c(0.3, 0.3)
fit <- fitCopulaOne(par, dat)
```