#' Density of the Smooth Additive-Perturbed Normal (SAPN) distribution
#'
#' @param x Numeric vector of quantiles.
#' @param mu Location parameter.
#' @param sigma Scale parameter, sigma > 0.
#' @param alpha Perturbation parameter, 0 < alpha <= 2.
#' @param c Perturbation location parameter.
#' @param k Steepness parameter of the sigmoid perturbation; default 2.
#' @return A numeric vector of density values, same length as `x`.
#' @note `mu`, `sigma`, `alpha`, `c`, and `k` must each be scalars (length 1).
#' @examples
#' dsapn(0, mu = 0, sigma = 1, alpha = 0.5, c = 2)
#' @export
dsapn <- function(x, mu, sigma, alpha, c, k = 2) {
  stopifnot(length(mu) == 1, length(sigma) == 1, length(alpha) == 1,
            length(c) == 1, length(k) == 1)
  stopifnot(sigma > 0, alpha > 0, alpha <= 2)
  K <- 2 / (2 + alpha)
  base_x    <- dnorm(x, mean = mu, sd = sigma)
  sigmoid_x <- 1 / (1 + exp(-k * ((x - c) / sigma)))
  bump_x    <- (alpha / sigma^2) * dnorm(x, mean = c, sd = sigma) * (x - c)^2 * sigmoid_x
  K * (base_x + bump_x)
}

#' Cumulative distribution function of the SAPN distribution
#'
#' Computed by numerical integration of \code{\link{dsapn}}.
#'
#' @param q Numeric vector of quantiles.
#' @param mu Location parameter.
#' @param sigma Scale parameter, sigma > 0.
#' @param alpha Perturbation parameter, 0 < alpha <= 2.
#' @param c Perturbation location parameter.
#' @param k Steepness parameter of the sigmoid perturbation; default 2.
#' @return A numeric vector of probabilities, same length as `q`.
#' @examples
#' psapn(1, mu = 0, sigma = 1, alpha = 0.5, c = 2)
#' @export
psapn <- function(q, mu, sigma, alpha, c, k = 2) {

  stopifnot(
    length(mu) == 1,
    length(sigma) == 1,
    length(alpha) == 1,
    length(c) == 1,
    length(k) == 1
  )

  stopifnot(sigma > 0, alpha > 0, alpha <= 2)

  vapply(q, function(v) {

    if (v >= max(mu, c) + 10 * sigma) {

      upper_tail <- tryCatch(
        integrate(
          dsapn,
          lower = v,
          upper = Inf,
          mu = mu,
          sigma = sigma,
          alpha = alpha,
          c = c,
          k = k
        )$value,
        error = function(e) NA_real_
      )

      return(1 - upper_tail)

    } else {

      tryCatch(
        integrate(
          dsapn,
          lower = -Inf,
          upper = v,
          mu = mu,
          sigma = sigma,
          alpha = alpha,
          c = c,
          k = k
        )$value,
        error = function(e) NA_real_
      )
    }

  }, numeric(1))
}

# Internal: builds a monotone-spline inverse-CDF via trapezoidal integration
# of dsapn on a wide grid. Shared engine for qsapn() and rsapn().
.sapn_build_sampler <- function(mu, sigma, alpha, c, k = 2,
                                grid_width = 20, n_grid = 5000) {

  lo <- min(mu, c) - grid_width * sigma
  hi <- max(mu, c) + grid_width * sigma
  x_grid <- seq(lo, hi, length.out = n_grid)
  f_vals <- dsapn(x_grid, mu, sigma, alpha, c, k)

  dx     <- diff(x_grid)
  area   <- 0.5 * (f_vals[-1] + f_vals[-n_grid]) * dx
  F_vals <- c(0, cumsum(area))
  F_vals <- F_vals / F_vals[n_grid]           # force F(hi) = 1 exactly

  keep   <- !duplicated(F_vals)               # guard flat tail regions
  F_vals <- F_vals[keep]
  x_grid <- x_grid[keep]

  splinefun(x = F_vals, y = x_grid, method = "monoH.FC")
}

#' Quantile function of the SAPN distribution
#'
#' @param p Numeric vector of probabilities.
#' @param mu Location parameter.
#' @param sigma Scale parameter, sigma > 0.
#' @param alpha Perturbation parameter, 0 < alpha <= 2.
#' @param c Perturbation location parameter.
#' @param k Steepness parameter of the sigmoid perturbation; default 2.
#' @param grid_width Half-width of the evaluation grid, in multiples of sigma; default 20.
#' @param n_grid Number of grid points used to build the inverse-CDF spline; default 5000.
#' @return A numeric vector of quantiles, same length as `p`.
#' @examples
#' qsapn(0.5, mu = 0, sigma = 1, alpha = 0.5, c = 2)
#' @export
qsapn <- function(p, mu, sigma, alpha, c, k = 2,
                  grid_width = 20, n_grid = 5000) {

  stopifnot(
    length(mu) == 1, length(sigma) == 1, length(alpha) == 1,
    length(c) == 1, length(k) == 1
  )
  stopifnot(sigma > 0, alpha > 0, alpha <= 2)

  # Setup output vector and identify valid probabilities
  out <- rep(NaN, length(p))
  valid <- !is.na(p) & p >= 0 & p <= 1

  if (any(!valid & !is.na(p))) {
    warning("NaNs produced: probabilities must lie in [0, 1].", call. = FALSE)
  }

  # Only build the spline and compute if there are valid probabilities
  if (any(valid)) {
    inv_cdf <- .sapn_build_sampler(
      mu, sigma, alpha, c, k, grid_width, n_grid
    )
    # Prevent edge extrapolation by clamping strictly inside (0, 1)
    p_clamped <- pmin(pmax(p[valid], 1e-10), 1 - 1e-10)
    out[valid] <- inv_cdf(p_clamped)
  }

  out
}

#' Random generation from the SAPN distribution
#'
#' @param n Number of observations to generate. If length(n) > 1, the length is used.
#' @param mu Location parameter.
#' @param sigma Scale parameter, sigma > 0.
#' @param alpha Perturbation parameter, 0 < alpha <= 2.
#' @param c Perturbation location parameter.
#' @param k Steepness parameter of the sigmoid perturbation; default 2.
#' @param grid_width Half-width of the evaluation grid, in multiples of sigma; default 20.
#' @param n_grid Number of grid points used to build the inverse-CDF spline; default 5000.
#' @return A numeric vector of length `n`.
#' @examples
#' rsapn(10, mu = 0, sigma = 1, alpha = 0.5, c = 2)
#' @export
rsapn <- function(n, mu, sigma, alpha, c, k = 2,
                  grid_width = 20, n_grid = 5000) {

  # Standard base-R handling of `n`
  if (length(n) > 1) n <- length(n)
  stopifnot(length(n) == 1, n >= 0)
  if (n == 0) return(numeric(0))

  stopifnot(length(mu) == 1, length(sigma) == 1, length(alpha) == 1,
            length(c) == 1, length(k) == 1)

  inv_cdf <- .sapn_build_sampler(mu, sigma, alpha, c, k, grid_width, n_grid)

  u <- runif(n)
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
  inv_cdf(u)
}
