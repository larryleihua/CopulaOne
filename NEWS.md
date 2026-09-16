# CopulaOne 0.0.1.201

- Exported and documented the FRA1 CDF, log-density, conditional CDF/inverse,
  generators, Kendall distribution, and tail summaries.
- Moved the main FRA1 calculations onto the log scale. Clipping is opt-in
  (`eps > 0`); default `eps = 0` preserves valid inputs and rejects invalid ones.
- Marginal quantiles have exact endpoints and checked convergence.
- Repaired likelihood validation, parameter grouping, bounded optimization,
  worker cleanup, and Hessian standard errors. Worker count defaults to 1.
- Dependence integrations report failures instead of substituting zeros.
- Repaired PPPP/GGEE extreme-value special cases and tied/missing uniform scores.
- Removed unbounded gamma resampling loops; GGGG remains internal.
- Seeded simulations may differ from earlier development versions because
  the random-generation algorithms now avoid intermediate underflow.
- Added regression tests, updated package documentation, and Windows CI checks.
