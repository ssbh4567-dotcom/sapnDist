# Internal: Anderson-Darling statistic of `data` against the SAPN CDF at given parameters
.sapn_ad_stat <- function(data, mu, sigma, alpha, c, k = 2) {
  as.numeric(
    ad.test(data, null = psapn,
            mu = mu, sigma = sigma, alpha = alpha, c = c, k = k,
            estimated = FALSE)$statistic
  )
}

#' Parametric bootstrap goodness-of-fit test for a fitted SAPN distribution
#'
#' @param x Original data vector.
#' @param fit An object of class `sapn_fit`, from \code{\link{fit_sapn_mme}}
#'   or \code{\link{fit_sapn_mle}}.
#' @param B Number of bootstrap replications; default 500.
#' @param seed Integer seed for reproducibility; default 123.
#' @param n_starts_boot Number of multi-start optimizer restarts used when
#'   refitting each bootstrap replicate; default 200.
#' @return A list with elements `statistic` (observed Anderson-Darling
#'   statistic), `p_value` (bootstrap p-value), `B_valid` (number of
#'   successful bootstrap replicates), and `method`.
#' @importFrom goftest ad.test
#' @examples
#' \donttest{
#' x <- rsapn(100, mu = 0, sigma = 1, alpha = 0.5, c = 2)
#' fit <- fit_sapn_mle(x, n_starts = 50)
#' gof_sapn(x, fit, B = 20, n_starts_boot = 20)
#' }
#' @export
gof_sapn <- function(x, fit, B = 500, seed = 123, n_starts_boot = 200) {
  if (!inherits(fit, "sapn_fit"))
    stop("`fit` must be a sapn_fit object from fit_sapn_mme() or fit_sapn_mle().")

  set.seed(seed)
  p  <- as.list(fit$estimate)
  A0 <- .sapn_ad_stat(x, p$mu, p$sigma, p$alpha, p$c, k = fit$k)

  refit <- function(xx) {
    if (identical(fit$method, "MLE")) {
      fit_sapn_mle(xx, k = fit$k, n_starts = n_starts_boot, hessian = FALSE)
    } else {
      fit_sapn_mme(xx, k = fit$k, n_starts = n_starts_boot)
    }
  }

  n <- length(x)
  A_boot <- numeric(B)
  for (b in seq_len(B)) {
    x_star   <- rsapn(n, p$mu, p$sigma, p$alpha, p$c, k = fit$k)
    fit_star <- try(refit(x_star), silent = TRUE)
    A_boot[b] <- if (inherits(fit_star, "try-error")) {
      NA_real_
    } else {
      ps <- as.list(fit_star$estimate)
      tryCatch(.sapn_ad_stat(x_star, ps$mu, ps$sigma, ps$alpha, ps$c, k = fit$k),
               error = function(e) NA_real_)
    }
  }
  A_boot <- A_boot[!is.na(A_boot)]

  list(
    statistic = A0,
    p_value   = (1 + sum(A_boot >= A0)) / (length(A_boot) + 1),
    B_valid   = length(A_boot),
    method    = fit$method
  )
}





