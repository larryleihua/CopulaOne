# CopulaOne - an R package for a very flexible bivariate parametric copula

The R package *CopulaOne* implements functions for a bivariate copula that is very flexible and parsimonious. 

Bivariate copulas have been widely used either in modeling bivariate dependence structures or building multivariate dependence models such as Vine copulas and factor copulas. In the literature, there are numurous parametric bivariate copula families. It is often very time consuming to select copula families from many different candidate copula families. The R package CopulaOne aims at implementing a very flexible bivariate copula that is parsimonious but very flexible and has the capicity of capturing most bivariate dependence patterns. Compared to those exisitng bivariate parametric copula families, the main merits of the bivariate copula is that, it can account for full-range tail dependence in both upper and lower tails, and the upper and lower tails can be either symmetric or asymmetric. 

## Installation
- The R package CopulaOne can be easily installed from github by the following two lines.
```{r, eval=FALSE}
library(devtools)
install_github("larryleihua/CopulaOne")
```

## Basic functions
- Simulation based on the copula can be done follows:
```{r}
library(CopulaOne)
UU = rGGEE_COP(10, a=0.5, b=0.8)
print(t(round(UU,4)))
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

## Model fitting
- An example of fitting dependence between exachange rates returns
```{r}
data("euro0306")
dat <- uscore(euro0306[,c(2,3)])[1:50,]
par <- c(0.3, 0.3)
fit <- fitCopulaOne(par, dat)
```