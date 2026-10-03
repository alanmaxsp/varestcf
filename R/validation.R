.need <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

.scalar_positive <- function(x, name) {
  .need(is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0,
        paste0(name, " must be a finite positive number."))
}

.ids <- function(x, n, name) {
  .need(is.character(x) && length(x) == n && !anyNA(x) &&
          all(nzchar(x)) && !anyDuplicated(x),
        paste0(name, " must contain unique, nonempty character identifiers."))
}

.matrix_input <- function(x, name) {
  .need((is.matrix(x) && is.numeric(x) && !is.complex(x)) ||
          methods::is(x, "dMatrix"),
        paste0(name, " must be a numeric matrix or a double Matrix object."))
  values <- if (methods::is(x, "dMatrix")) x@x else x
  .need(all(is.finite(values)), paste0(name, " contains missing or nonfinite values."))
}

.sparse <- function(x) {
  methods::as(methods::as(Matrix::Matrix(x, sparse = TRUE), "generalMatrix"),
              "CsparseMatrix")
}

.named_square <- function(x, ids, name) {
  .matrix_input(x, name)
  .need(nrow(x) == length(ids) && ncol(x) == length(ids),
        paste0(name, " has incompatible dimensions."))
  .ids(rownames(x), length(ids), paste0(name, " row names"))
  .ids(colnames(x), length(ids), paste0(name, " column names"))
  .need(setequal(rownames(x), ids) && setequal(colnames(x), ids),
        paste0(name, " identifiers do not match the required identifiers."))
  x[ids, ids, drop = FALSE]
}

.symmetric <- function(x, name, tolerance) {
  scale <- max(abs(x))
  tx <- if (methods::is(x, "Matrix")) Matrix::t(x) else t(x)
  asymmetry <- max(abs(x - tx))
  .need(asymmetry <= tolerance * max(scale, .Machine$double.xmin),
        paste0(name, " is not symmetric within tolerance."))
  # Only average discrepancies already established to be numerical roundoff.
  (x + tx) / 2
}

.spd_factor <- function(x, name) {
  result <- tryCatch(
    withCallingHandlers(
      Matrix::Cholesky(Matrix::forceSymmetric(x, uplo = "U"),
                       LDL = FALSE, perm = TRUE, super = TRUE),
      warning = function(w) stop(conditionMessage(w), call. = FALSE)),
    error = function(e) stop(paste0(name,
      " could not be factored as positive definite. Check identification, ",
      "precision matrices and numerical conditioning. Details: ", conditionMessage(e)),
      call. = FALSE))
  result
}

.memory_guard <- function(bytes, budget, stage) {
  .need(is.finite(bytes) && bytes <= budget * 1024^3,
        paste0(stage, " exceeds memory_budget_gib (estimated working storage ",
               format(round(bytes / 1024^3, 4), trim = TRUE),
               " GiB). This budget is an advisory guard, not an OS memory limit."))
}
