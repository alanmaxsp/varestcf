.varest_from_covariance <- function(ebv, Sigma, normalization = c("n", "n-1"), tolerance = 1e-8) {
  started <- proc.time()[[3L]]
  normalization <- match.arg(normalization)
  .scalar_positive(tolerance, "tolerance")
  .need(is.numeric(ebv) && is.null(dim(ebv)) && !is.complex(ebv) &&
          length(ebv) >= 2L && all(is.finite(ebv)),
        "ebv must be a finite numeric vector containing at least two values.")
  n <- length(ebv); ids <- names(ebv)
  .ids(ids, n, "ebv names")
  .need(!methods::is(Sigma, "sparseMatrix"),
        "Sigma must be an explicitly complete dense matrix; sparse omissions cannot certify coverage.")
  S <- .named_square(Sigma, ids, "Sigma")
  # Full focal matrices are intentional. Missing entries must be NA, never implicit zeros.
  S <- as.matrix(S)
  S <- .symmetric(S, "Sigma", tolerance)
  scale <- max(abs(S))
  positive <- tryCatch({chol(S); TRUE}, error = function(e) FALSE)
  min_eigenvalue <- NA_real_
  if (!positive) {
    # Semidefinite external covariances are valid (for example, a zero error matrix).
    min_eigenvalue <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
    .need(min_eigenvalue >= -tolerance * max(scale, .Machine$double.xmin),
          "Sigma is not positive semidefinite within tolerance.")
  }
  centered <- ebv - mean(ebv)
  quadratic <- sum(centered^2)
  raw_trace <- sum(diag(S)) - sum(S) / n
  trace_scale <- sum(abs(diag(S))) + sum(abs(S)) / n
  .need(raw_trace >= -tolerance * max(trace_scale, .Machine$double.xmin),
        "Sigma has a negative centered error correction beyond tolerance.")
  centered_trace <- max(0, raw_trace)
  denominator <- if (normalization == "n") n else n - 1L
  estimate <- (quadratic + centered_trace) / denominator
  .need(is.finite(estimate), "VarEst-CF overflowed; check numerical scale.")
  list(estimate = estimate, ebv_term = quadratic / denominator,
       error_term = centered_trace / denominator, n = n,
       normalization = normalization, denominator = denominator,
       diagnostics = list(ids = ids, raw_centered_trace = raw_trace,
         trace_roundoff_clamped = raw_trace < 0, covariance_check = if (positive)
           "positive definite (Cholesky)" else "positive semidefinite (eigenvalues)",
         minimum_eigenvalue_if_computed = min_eigenvalue, tolerance = tolerance,
         seconds = proc.time()[[3L]] - started))
}
