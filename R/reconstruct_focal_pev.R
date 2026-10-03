# Reference implementation: C contains the full identified fixed/random model.
# Q is the precision of the RANDOM EFFECTS, including their variance scale.
.reconstruct_from_matrices <- function(X, Z, random_precision, focal_ids,
                                  residual_precision = 1, block_size = 8L,
                                  memory_budget_gib = 1, tolerance = 1e-8, y = NULL) {
  started <- proc.time()[[3L]]
  .scalar_positive(tolerance, "tolerance")
  .scalar_positive(memory_budget_gib, "memory_budget_gib")
  .scalar_positive(block_size, "block_size")
  .need(block_size == floor(block_size) && block_size <= .Machine$integer.max,
        "block_size must be a positive integer within the R integer range.")
  .matrix_input(X, "X"); .matrix_input(Z, "Z")
  n <- nrow(Z); p <- ncol(X); r <- ncol(Z)
  .need(n > 0L && r > 0L && nrow(X) == n,
        "X and Z must have the same positive number of records; Z needs columns.")
  .ids(colnames(Z), r, "Z column names")
  ids <- colnames(Z)
  .ids(focal_ids, length(focal_ids), "focal_ids")
  .need(length(focal_ids) > 0L && all(focal_ids %in% ids),
        "focal_ids must select at least one coefficient present in Z.")
  # Record ordering is either explicitly named in both designs or positional in both.
  if (!is.null(rownames(X)) || !is.null(rownames(Z))) {
    .ids(rownames(X), n, "X row names"); .ids(rownames(Z), n, "Z row names")
    .need(identical(rownames(X), rownames(Z)), "X and Z record orders differ.")
  }
  q <- length(focal_ids); m <- as.double(p) + r
  .need(m < .Machine$integer.max, "The equation dimension exceeds the CSC index limit.")
  block_size <- as.integer(min(block_size, q))
  # Output, symmetry/validation work and solve/residual buffers; the factor is unknown.
  output_work <- 32 * as.double(q)^2
  buffers <- 40 * m * block_size
  .memory_guard(output_work + buffers, memory_budget_gib, "Focal output and buffers")

  Q <- .named_square(random_precision, ids, "random_precision")
  X <- .sparse(X); Z <- .sparse(Z)
  Q <- .sparse(.symmetric(.sparse(Q), "random_precision", tolerance))
  W <- residual_precision
  if (is.numeric(W) && is.null(dim(W))) {
    .need(length(W) %in% c(1L, n) && all(is.finite(W)) && all(W > 0),
          "residual_precision must be positive and scalar or one value per record.")
    if (!is.null(names(W))) {
      .need(length(W) == n && !is.null(rownames(Z)),
            "Named residual precisions require named records and one value per record.")
      .ids(names(W), n, "residual_precision names")
      .need(setequal(names(W), rownames(Z)), "Residual precision record identifiers differ.")
      W <- W[rownames(Z)]
    }
    W <- Matrix::Diagonal(n, x = rep_len(as.numeric(W), n))
    residual_kind <- "diagonal"
  } else {
    .need(!is.null(rownames(Z)), "A residual precision matrix requires named records.")
    W <- .named_square(W, rownames(Z), "residual_precision")
    W <- .sparse(.symmetric(.sparse(W), "residual_precision", tolerance))
    residual_kind <- "matrix"
  }
  input_bytes <- sum(vapply(list(X, Z, Q, W), function(a) as.numeric(object.size(a)), 0))
  .memory_guard(input_bytes + output_work + buffers, memory_budget_gib,
                "Inputs, focal output and buffers")
  # Strict reference checks: this incurs additional factorizations, reported separately.
  check_factor <- .spd_factor(Q, "random_precision")
  .memory_guard(input_bytes + as.numeric(object.size(check_factor)) + output_work + buffers,
                memory_budget_gib, "Random precision validation")
  rm(check_factor)
  if (residual_kind == "matrix") {
    check_factor <- .spd_factor(W, "residual_precision")
    .memory_guard(input_bytes + as.numeric(object.size(check_factor)) + output_work + buffers,
                  memory_budget_gib, "Residual precision validation")
    rm(check_factor)
  }
  checked <- proc.time()[[3L]]
  D <- cbind(X, Z)
  rhs <- NULL
  if (!is.null(y)) {
    .need(is.numeric(y) && is.null(dim(y)) && length(y) == n && all(is.finite(y)),
          "The response must contain one finite numeric value per model record.")
    rhs <- as.numeric(Matrix::crossprod(D, W %*% y))
  }
  # General CSC products cannot store >= INT_MAX positions. This structural bound
  # can overestimate storage but prevents knowingly constructing oversized products.
  if (residual_kind == "diagonal") {
    degree <- as.numeric(Matrix::rowSums(D != 0))
    data_bound <- min(m^2, sum(degree^2))
  } else {
    # Correlated residuals can connect every pair of coefficients.
    data_bound <- m^2
  }
  union_bound <- min(m^2, data_bound + length(Q@x))
  .need(union_bound < .Machine$integer.max,
        "The conservative assembly bound exceeds the 32-bit CSC index limit.")
  assembly_bound <- input_bytes + 24 * union_bound + 8 * (m + 1) + output_work + buffers
  .memory_guard(assembly_bound, memory_budget_gib, "Conservative assembly bound")
  penalty <- if (p) Matrix::bdiag(Matrix::Diagonal(p, x = 0), Q) else Q
  C <- Matrix::forceSymmetric(Matrix::crossprod(D, W %*% D) + penalty, uplo = "U")
  assembled <- proc.time()[[3L]]
  C_bytes <- as.numeric(object.size(C))
  rm(D, penalty, X, Z, Q, W)
  factor <- .spd_factor(C, "The mixed-model equations (X must have full column rank)")
  factored <- proc.time()[[3L]]
  factor_bytes <- as.numeric(object.size(factor))
  .memory_guard(C_bytes + factor_bytes + output_work + buffers, memory_budget_gib,
                "Factored equations and planned extraction")
  indices <- p + match(focal_ids, ids)
  ebv <- NULL
  if (!is.null(rhs)) {
    solution <- as.numeric(Matrix::solve(factor, rhs, system = "A"))
    error <- max(abs(as.numeric(C %*% solution) - rhs))
    denominator <- as.numeric(Matrix::norm(C, "I")) * max(abs(solution)) + max(abs(rhs))
    .need(is.finite(error) && is.finite(denominator) &&
            error <= tolerance * max(denominator, .Machine$double.xmin),
          "The EBV solution failed the backward-error check.")
    ebv <- setNames(solution[indices], focal_ids)
  }
  Sigma <- matrix(NA_real_, q, q, dimnames = list(focal_ids, focal_ids))
  starts <- seq.int(1L, q, by = block_size)
  checks <- vector("list", length(starts))
  C_norm <- as.numeric(Matrix::norm(C, type = "I"))
  for (b in seq_along(starts)) {
    columns <- seq.int(starts[b], min(q, starts[b] + block_size - 1L))
    E <- matrix(0, m, length(columns))
    E[cbind(indices[columns], seq_along(columns))] <- 1
    V <- as.matrix(Matrix::solve(factor, E, system = "A"))
    residual <- max(abs(as.matrix(C %*% V) - E))
    denominator <- C_norm * max(abs(V)) + 1
    .need(is.finite(denominator) && is.finite(residual), "Nonfinite solve diagnostics.")
    backward_error <- residual / denominator
    .need(backward_error <= tolerance, "Focal solve failed the backward-error tolerance.")
    Sigma[, columns] <- V[indices, , drop = FALSE]
    checks[[b]] <- data.frame(first = min(columns), last = max(columns),
                             absolute_residual = residual, backward_error = backward_error)
  }
  solved <- proc.time()[[3L]]
  .need(all(is.finite(Sigma)), "Focal covariance contains nonfinite values.")
  asymmetry <- max(abs(Sigma - t(Sigma)))
  Sigma <- .symmetric(Sigma, "Reconstructed focal covariance", tolerance)
  tryCatch(chol(Sigma), error = function(e)
    stop("Reconstructed focal covariance is not numerically positive definite.", call. = FALSE))
  finished <- proc.time()[[3L]]
  list(Sigma = Sigma, focal_ids = focal_ids, ebv = ebv,
       diagnostics = list(records = n, fixed_coefficients = p, random_coefficients = r,
         equations = m, focal_coefficients = q, block_size = block_size,
         precision_checks = "positive definiteness checked by factorization",
         C_stored_entries = length(C@x), C_bytes = C_bytes, factor_bytes = factor_bytes,
         Sigma_bytes = as.numeric(object.size(Sigma)), output_values_bytes = 8 * as.double(q)^2,
         assembly_bound_bytes = assembly_bound, memory_budget_gib = memory_budget_gib,
         memory_note = "Advisory estimates/object sizes; not measured peak RSS or an OS limit.",
         max_asymmetry = asymmetry, tolerance = tolerance, solves = do.call(rbind, checks),
         seconds = c(validation = checked - started, assembly = assembled - checked,
                     factorization = factored - assembled, extraction = solved - factored,
                     final_validation = finished - solved, total = finished - started)))
}
