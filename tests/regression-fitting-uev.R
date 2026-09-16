pkg <- asNamespace("CopulaOne")
worker_count <- get(".copula_worker_count", pkg)
complete_nllk <- get(".complete_copula_nllk", pkg)
seq_run <- get("seqRun", pkg)
pppp_density <- get("dPPPP_COP", pkg)
uev_cdf <- get("pUEV_GGEE_COP", pkg)
uev_density <- get("dUEV_GGEE_COP", pkg)

stopifnot(worker_count(2, 64) == 2L,
          worker_count(2, 1) == 1L,
          worker_count(2, NA_integer_) == 1L)

# Every observation must contribute exactly once, even with empty chunks.
dat <- cbind(c(.2, .35, .7), c(.3, .6, .8))
expected <- -log(pppp_density(dat[, 1], dat[, 2], 1, 1, 1, 1))
for (workers in c(1L, 2L, 3L, 5L)) {
  parts <- lapply(seq_len(workers), seq_run, dat = dat, nco = workers,
                  para = rep(1, 4), copula_family = "PPPP")
  stopifnot(isTRUE(all.equal(unlist(parts), expected)),
            isTRUE(all.equal(complete_nllk(parts, nrow(dat)), sum(expected))))
}
single <- seq_run(1L, dat[1, , drop = FALSE], 1L, rep(1, 4))
stopifnot(isTRUE(all.equal(single, expected[1])))

# Invalid, missing, and failed likelihood evaluations must be penalized.
for (parts in list(list(NA_real_, 2), list(Inf, 2), list(NaN, 2),
                   list(-Inf, 2), list(2), list(),
                   structure("worker failed", class = "try-error"))) {
  stopifnot(complete_nllk(parts, 2L) == 1e100)
}
stopifnot(complete_nllk(list(-2, 3), 2L) == 1,
          complete_nllk(list(1e308, 1e308), 2L) == 1e100)

# Outside the open parameter interval, both functions describe independence.
for (b in c(-Inf, -1, 0, 1, 1.2, Inf)) {
  stopifnot(uev_cdf(.9, .9, b) == .81,
            uev_cdf(.3, .4, b) == .12,
            uev_density(.3, .4, b) == 1,
            uev_density(.9, .9, b) == 1)
}
for (b in c(.1, .5, .9)) {
  stopifnot(uev_cdf(0, .4, b) == 0,
            uev_cdf(.4, 0, b) == 0,
            uev_cdf(1, .4, b) == .4,
            uev_cdf(.4, 1, b) == .4,
            uev_cdf(.9, .9, b) == .9^(1 + b))
  for (u in c(.1, .3, .7, .9)) for (v in c(.1, .3, .7, .9)) {
    value <- uev_cdf(u, v, b)
    stopifnot(is.finite(value), value >= max(u + v - 1, 0) - 1e-12,
              value <= min(u, v) + 1e-12)
  }
}
