pkg <- asNamespace("CopulaOne")
cdf <- get("pUEV_PPPP_COP", pkg)
density <- get("dUEV_PPPP_COP", pkg)
score <- get("uscore", pkg)
near <- function(x, y, tol = 1e-10)
  stopifnot(isTRUE(all.equal(x, y, tolerance = tol)))

# An independent expression using the corrected b^2 coefficient, away
# from the removable singularity. This catches a regression to a^2.
reference_cdf <- function(u, v, be, a, b) {
  x <- min(-log(u), -log(v)); y <- max(-log(u), -log(v))
  A1 <- be * (a/(a+b)/(be-b) - a*b/(a+b)^2/(a-b+be) +
                (a/(a+b))^2/(2*b-be))
  A2 <- be * (a*b/(a+b)^2/(a-b+be) - b/(a+b)/(a+be) +
                (b/(a+b))^2/(2*a+be))
  A3 <- be * (1/be - b/(a+b)/(a+be) - a/(a+b)/(be-b))
  exp(-x-y + (a+be)*(b-be)/(a*b) *
        (A1*x^(b/be)*y^(1-b/be) + A2*x^(1+a/be)*y^(-a/be) + A3*x))
}
for (pars in list(c(.6, 2, 1), c(.4, .3, 2), c(.6, 1, 1))) {
  near(cdf(.3, .4, pars[1], pars[2], pars[3]),
       reference_cdf(.3, .4, pars[1], pars[2], pars[3]))
}

for (be in c(1, 1.2, 2)) {
  near(cdf(.3, .4, be, 1, 1), .12)
  near(density(.3, .4, be, 1, 1), 1)
  near(density(.5, .5, be, 1, 1), 1)
}

# Check the singular limit, both sides of it, and the independence limit.
for (fun in list(cdf, density)) {
  exact <- fun(.3, .4, 1, 1, 2)
  stopifnot(is.finite(exact))
  for (offset in c(-1e-8, 1e-8))
    near(fun(.3, .4, 1, 1, 2 + offset), exact, 1e-7)
}
near(cdf(.3, .4, 1, 1, 1 + 1e-9), .12, 1e-7)
near(density(.3, .4, 1, 1, 1 + 1e-9), 1, 1e-7)

for (pars in list(c(.6, 2, 1), c(1, 1, 2), c(.4, .3, 2))) {
  C <- function(u, v) cdf(u, v, pars[1], pars[2], pars[3])
  d <- function(u, v) density(u, v, pars[1], pars[2], pars[3])
  near(C(0, .4), 0); near(C(.4, 0), 0)
  near(C(1, .4), .4); near(C(.4, 1), .4)
  for (u in c(.1, .3, .7, .9)) for (v in c(.1, .3, .7, .9)) {
    value <- C(u, v)
    stopifnot(is.finite(value), value >= max(u + v - 1, 0) - 1e-12,
              value <= min(u, v) + 1e-12, is.finite(d(u, v)), d(u, v) >= 0)
    near(value, C(v, u)); near(d(u, v), d(v, u))
    h <- 1e-4
    mixed <- (C(u+h, v+h) - C(u+h, v-h) - C(u-h, v+h) + C(u-h, v-h))/(4*h^2)
    near(d(u, v), mixed, 1e-4)
  }
  # A copula density must integrate to one along either margin.
  for (u in c(.2, .5, .8)) {
    integral <- stats::integrate(function(v) vapply(v, function(vv) d(u, vv),
                                                   numeric(1)), 0, 1,
                                 rel.tol = 1e-6, subdivisions = 500L)
    near(integral$value, 1, 1e-5)
  }
}

near(score(c(1, 1, 2, 3)), c(.25, .25, .625, .875))
near(score(c(1, NA, 3)), c(.25, NA, .75))
near(score(c(3, 1, 2)), c(5/6, 1/6, .5))
near(score(c(1, 1, 2), aunif = 0), c(.375, .375, .75))
near(score(c(NA, NA)), c(NA_real_, NA_real_))
near(score(numeric(0)), numeric(0))
near(score(c(NA, 5, NA)), c(NA, .5, NA))
near(score(c(a = 1, b = 1)), c(a = .5, b = .5))
input <- data.frame(x = c(1, 1, NA, 3), y = c(NA, 2, 4, NA))
expected <- cbind(x = c(1/3, 1/3, NA, 5/6), y = c(NA, .25, .75, NA))
rownames(expected) <- rownames(input)
near(score(input), expected)
dimnames(expected) <- dimnames(as.matrix(input))
near(score(as.matrix(input)), expected)
stopifnot(identical(dim(score(matrix(numeric(0), 0, 2))), c(0L, 2L)))
