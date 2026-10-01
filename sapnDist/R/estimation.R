# Internal: negative log-likelihood
.sapn_nll <- function(params, data, k = 2) {
  mu <- params[1]; sigma <- params[2]; alpha <- params[3]; c <- params[4]
  pdf_vals <- dsapn(data, mu, sigma, alpha, c, k = k)
  pdf_vals <- pmax(pdf_vals, 1e-30)
  -sum(log(pdf_vals))
}

# Internal: analytical gradient of the NLL (exact score, used by L-BFGS-B)
.sapn_grad <- function(params, data, k = 2) {
  mu <- params[1]; sigma <- params[2]; alpha <- params[3]; c <- params[4]
  n  <- length(data)

  u_i   <- dnorm(data, mean = mu, sd = sigma)
  phi_c <- dnorm(data, mean = c,  sd = sigma)
  S_i   <- 1 / (1 + exp(-k * ((data - c) / sigma)))
  v_i   <- phi_c * (data - c)^2 * S_i

  denom <- u_i + (alpha / sigma^2) * v_i
  denom <- pmax(denom, 1e-30)

  d_mu <- (1 / sigma^2) * sum((u_i * (data - mu)) / denom)

  term_sigma <- (data - c)^2 - 2 * sigma^2 - k * sigma * (data - c) * (1 - S_i)
  d_sigma <- -(n / sigma) +
    (1 / sigma^3) * sum((u_i * (data - mu)^2 + (alpha / sigma^2) * v_i * term_sigma) / denom)

  d_alpha <- -(n / (2 + alpha)) + (1 / sigma^2) * sum(v_i / denom)

  term_c <- ((data - c)^2 / sigma^2) - 2 - (k / sigma) * (data - c) * (1 - S_i)
  d_c <- (alpha / sigma^2) * sum((phi_c * (data - c) * S_i * term_c) / denom)

  -c(d_mu, d_sigma, d_alpha, d_c)
}

# Internal: theoretical raw moments E[X^r], r = 1..4, via the Omega() integrals
.sapn_moments <- function(mu, sigma, alpha, c, k = 2) {
  K <- 2 / (2 + alpha)
  N1 <- mu; N2 <- mu^2 + sigma^2
  N3 <- mu^3 + 3 * mu * sigma^2
  N4 <- mu^4 + 6 * mu^2 * sigma^2 + 3 * sigma^4
  M2 <- sigma^2; M4 <- 3 * sigma^4; M6 <- 15 * sigma^6

  Omega <- function(m) {
    integrand <- function(x) (x^m) * dnorm(x, mean = 0, sd = sigma) * tanh((k * x) / (2 * sigma))
    tryCatch(
      integrate(integrand, lower = 0, upper = 20 * sigma,
                rel.tol = 1e-8, abs.tol = 1e-8)$value,
      error = function(e) NA_real_
    )
  }
  O3 <- Omega(3); O5 <- Omega(5)
  if (is.na(O3) || is.na(O5)) return(list(E1 = NA, E2 = NA, E3 = NA, E4 = NA))

  B1 <- O3 + c * (M2 / 2)
  B2 <- (M4 / 2) + 2 * c * O3 + c^2 * (M2 / 2)
  B3 <- O5 + 3 * c * (M4 / 2) + 3 * c^2 * O3 + c^3 * (M2 / 2)
  B4 <- (M6 / 2) + 4 * c * O5 + 6 * c^2 * (M4 / 2) + 4 * c^3 * O3 + c^4 * (M2 / 2)

  list(
    E1 = K * (N1 + (alpha / sigma^2) * B1),
    E2 = K * (N2 + (alpha / sigma^2) * B2),
    E3 = K * (N3 + (alpha / sigma^2) * B3),
    E4 = K * (N4 + (alpha / sigma^2) * B4)
  )
}

#' Theoretical skewness and excess kurtosis implied by SAPN parameters
#'
#' Excess kurtosis (beta2 - 3) is reported, matching the convention used
#' elsewhere in the paper's descriptive-statistics table.
#'
#' @param mu Location parameter.
#' @param sigma Scale parameter, sigma > 0.
#' @param alpha Perturbation parameter, 0 < alpha <= 2.
#' @param c Perturbation location parameter.
#' @param k Steepness parameter of the sigmoid perturbation; default 2.
#' @return A named numeric vector with elements `skewness` and `excess_kurtosis`.
#' @examples
#' sapn_shape(mu = 0, sigma = 1, alpha = 0.5, c = 2)
#' @export
sapn_shape <- function(mu, sigma, alpha, c, k = 2) {
  m <- .sapn_moments(mu, sigma, alpha, c, k)
  if (anyNA(unlist(m))) return(c(skewness = NA_real_, excess_kurtosis = NA_real_))
  mu2 <- m$E2 - m$E1^2
  mu3 <- m$E3 - 3 * m$E1 * m$E2 + 2 * m$E1^3
  mu4 <- m$E4 - 4 * m$E1 * m$E3 + 6 * m$E1^2 * m$E2 - 3 * m$E1^4
  c(skewness = mu3 / mu2^1.5, excess_kurtosis = mu4 / mu2^2 - 3)
}

# Internal: MME loss (sum of squared standardized raw-moment errors)
.sapn_mme_loss <- function(params, data, k = 2) {
  mu <- params[1]; sigma <- params[2]; alpha <- params[3]; c <- params[4]
  if (sigma <= 0 || alpha <= 0 || alpha > 2) return(1e9)
  theo <- .sapn_moments(mu, sigma, alpha, c, k = k)
  if (is.na(theo$E1)) return(1e9)

  sd_data <- sd(data)
  m1 <- mean(data); m2 <- mean(data^2); m3 <- mean(data^3); m4 <- mean(data^4)

  ((theo$E1 - m1) / sd_data)^2 + ((theo$E2 - m2) / sd_data^2)^2 +
    ((theo$E3 - m3) / sd_data^3)^2 + ((theo$E4 - m4) / sd_data^4)^2
}

# Internal: shared multi-start grid for both estimators
.sapn_start_grid <- function(x, n_starts) {
  sd_x <- sd(x)
  data.frame(
    mu    = rnorm(n_starts, mean(x), 2 * sd_x),
    sigma = exp(runif(n_starts, log(0.02 * sd_x), log(20 * sd_x))),
    alpha = runif(n_starts, 0.0001, 2),
    c     = runif(n_starts, min(x) - 2 * sd_x, max(x) + 2 * sd_x)
  )
}

#' Method-of-moments estimation for the SAPN distribution
#'
#' @param x Numeric data vector.
#' @param k Steepness parameter of the sigmoid perturbation; default 2.
#' @param n_starts Number of random starting values for the multi-start optimizer; default 1000.
#' @param seed Optional integer seed for reproducibility of the starting grid.
#' @return An object of class `sapn_fit` with elements `estimate`, `loss`, `k`, `method`, `n`.
#' @examples
#' \donttest{
#' fit_sapn_mme(rsapn(100, mu = 0, sigma = 1, alpha = 0.5, c = 2), n_starts = 50)
#' }
#' @export
fit_sapn_mme <- function(x, k = 2, n_starts = 1000, seed = NULL) {
  stopifnot(length(k) == 1, length(n_starts) == 1)
  if (!is.null(seed)) set.seed(seed)
  start_grid <- .sapn_start_grid(x, n_starts)

  best_loss <- Inf; best_par <- rep(NA_real_, 4)
  for (i in seq_len(n_starts)) {
    fit <- try(
      optim(as.numeric(start_grid[i, ]), .sapn_mme_loss, data = x, k = k,
            method = "Nelder-Mead", control = list(maxit = 20000)),
      silent = TRUE
    )
    if (!inherits(fit, "try-error") && fit$convergence == 0 && fit$value < best_loss) {
      best_loss <- fit$value; best_par <- fit$par
    }
  }
  if (!is.finite(best_loss)) stop("MME optimization failed to converge from any starting value.")

  out <- list(estimate = setNames(best_par, c("mu", "sigma", "alpha", "c")),
              loss = best_loss, k = k, method = "MME", n = length(x))
  class(out) <- "sapn_fit"
  out
}

#' Maximum-likelihood estimation for the SAPN distribution
#'
#' @param x Numeric data vector.
#' @param k Steepness parameter of the sigmoid perturbation; default 2.
#' @param n_starts Number of random starting values for the multi-start optimizer; default 1000.
#' @param seed Optional integer seed for reproducibility of the starting grid.
#' @param hessian Logical; if TRUE, also fit once more at the optimum to obtain
#'   the Hessian for standard errors. Default TRUE.
#' @return An object of class `sapn_fit` with elements `estimate`, `loglik`, `k`,
#'   `method`, `n`, and (if `hessian = TRUE`) `hessian`.
#' @examples
#' \donttest{
#' fit_sapn_mle(rsapn(100, mu = 0, sigma = 1, alpha = 0.5, c = 2), n_starts = 50)
#' }
#' @export
fit_sapn_mle <- function(x, k = 2, n_starts = 1000, seed = NULL, hessian = TRUE) {
  stopifnot(length(k) == 1, length(n_starts) == 1)
  if (!is.null(seed)) set.seed(seed)
  start_grid <- .sapn_start_grid(x, n_starts)

  best_nll <- Inf; best_par <- rep(NA_real_, 4)
  for (i in seq_len(n_starts)) {
    fit <- try(
      optim(as.numeric(start_grid[i, ]), .sapn_nll, gr = .sapn_grad, data = x, k = k,
            method = "L-BFGS-B",
            lower = c(-Inf, 1e-6, 1e-6, -Inf), upper = c(Inf, Inf, 2, Inf)),
      silent = TRUE
    )
    if (!inherits(fit, "try-error") && fit$convergence == 0 && fit$value < best_nll) {
      best_nll <- fit$value; best_par <- fit$par
    }
  }
  if (!is.finite(best_nll)) stop("MLE optimization failed to converge from any starting value.")

  out <- list(estimate = setNames(best_par, c("mu", "sigma", "alpha", "c")),
              loglik = -best_nll, k = k, method = "MLE", n = length(x))

  if (hessian) {
    hfit <- optim(best_par, .sapn_nll, gr = .sapn_grad, data = x, k = k,
                  method = "L-BFGS-B",
                  lower = c(-Inf, 1e-6, 1e-6, -Inf), upper = c(Inf, Inf, 2, Inf),
                  hessian = TRUE)
    out$hessian <- hfit$hessian
  }
  class(out) <- "sapn_fit"
  out
}

#' Confidence intervals for a fitted SAPN model
#'
#' Wald and log-scale (delta-method) confidence intervals, computed from the
#' Hessian of a fit produced by \code{\link{fit_sapn_mle}}.
#'
#' @param object An object of class `sapn_fit`, from \code{\link{fit_sapn_mle}}
#'   with `hessian = TRUE`.
#' @param parm Optional character vector of parameter names ("mu", "sigma",
#'   "alpha", "c") to include. If missing, all four are returned.
#' @param level Confidence level; default 0.95.
#' @param ... Additional arguments, currently unused.
#' @return A data frame with one row per parameter, giving the estimate,
#'   standard error, and Wald and log-scale confidence bounds.
#' @note When the fitted `alpha` is close to the boundary `alpha = 0`
#'   (below \eqn{10^{-4}}), `c` becomes weakly identifiable and the observed
#'   Fisher information matrix may be nearly singular; a warning is issued
#'   in this case, and the resulting intervals should be interpreted with
#'   caution.
#' @importFrom stats confint qnorm
#' @export
confint.sapn_fit <- function(object, parm, level = 0.95, ...) {
  if (!identical(object$method, "MLE") || is.null(object$hessian))
    stop("confint.sapn_fit() requires a fit from fit_sapn_mle(..., hessian = TRUE).")

  alpha_hat <- unname(object$estimate["alpha"])
  if (alpha_hat < 1e-4) {
    warning(
      "alpha_hat = ", format(alpha_hat, scientific = TRUE),
      " is close to the boundary alpha = 0. Near this boundary, c becomes ",
      "weakly identifiable and the observed Fisher information matrix may ",
      "be nearly singular; the intervals below should be interpreted with ",
      "caution (see the Remark on Identifiability).",
      call. = FALSE
    )
  }

  cov_mat <- tryCatch(
    solve(object$hessian),
    error = function(e) {
      stop(
        "Hessian is computationally singular -- consistent with weak ",
        "identifiability near alpha = 0 (see the Remark on Identifiability). ",
        "Wald/delta-method confidence intervals are not available for this fit.",
        call. = FALSE
      )
    }
  )
  se  <- sqrt(diag(cov_mat))
  est <- object$estimate
  z   <- qnorm(1 - (1 - level) / 2)

  wald_lo <- est - z * se
  wald_hi <- est + z * se

  log_eligible <- names(est) %in% c("sigma", "alpha")
  log_lo <- ifelse(log_eligible, est * exp(-z * se / est), NA_real_)
  log_hi <- ifelse(log_eligible, est * exp( z * se / est), NA_real_)

  out <- data.frame(
    parameter  = names(est), estimate = unname(est), se = unname(se),
    wald_lower = unname(wald_lo),  wald_upper = unname(wald_hi),
    log_lower  = unname(log_lo),   log_upper  = unname(log_hi),
    row.names  = NULL
  )
  if (!missing(parm)) out <- out[out$parameter %in% parm, , drop = FALSE]
  out
}
#' Print a fitted SAPN model
#'
#' @param x An object of class `sapn_fit`.
#' @param ... Additional arguments, currently unused.
#' @return Invisibly returns `x`.
#' @export
print.sapn_fit <- function(x, ...) {
  cat("\nFitted SAPN Distribution (Method:", x$method, ")\n")
  cat("Observations:", x$n, "\n\n")
  cat("Parameter Estimates:\n")
  print(x$estimate)
  if (x$method == "MLE") {
    cat("\nLog-Likelihood:", x$loglik, "\n")
  } else {
    cat("\nMME Loss:", x$loss, "\n")
  }
  invisible(x)
}

