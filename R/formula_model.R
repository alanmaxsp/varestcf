.identified_design <- function(X, tolerance) {
  if (!ncol(X)) return(list(X = X, removed = character()))
  X <- .sparse(X)
  .need(all(is.finite(X@x)), "Fixed-effect predictors contain missing or nonfinite values.")
  scale <- sqrt(Matrix::colSums(X^2))
  active <- which(scale > 0)
  if (!length(active)) return(list(X = X[, FALSE, drop = FALSE], removed = colnames(X)))
  normalized <- X[, active, drop = FALSE] %*% Matrix::Diagonal(x = 1 / scale[active])
  if (nrow(normalized) < ncol(normalized))
    normalized <- rbind(normalized, Matrix::Matrix(0, ncol(normalized) - nrow(normalized), ncol(normalized), sparse = TRUE))
  qr <- Matrix::qr(normalized)
  keep <- active[qr@q[which(abs(Matrix::diag(qr@R)) > tolerance)] + 1L]
  keep <- sort(keep)
  list(X = X[, keep, drop = FALSE], removed = colnames(X)[setdiff(seq_len(ncol(X)), keep)])
}

.formula_model <- function(data, formula, animal, focal_ids,
                           relationship, relationship_inverse, random,
                           genetic_variance, residual_variance,
                           budget, tolerance, allow_missing = FALSE,
                           template = NULL, prepare_precision = TRUE) {
  .need(is.data.frame(data) && nrow(data) > 0L, "data must be a nonempty data.frame.")
  data <- as.data.frame(data)
  .need(inherits(formula, "formula") && length(formula) == 3L,
        "formula must include a response, for example weight ~ sex + group.")
  .need(is.character(animal) && length(animal) == 1L && animal %in% names(data),
        "animal must name the animal-ID column in data.")
  .need(inherits(random, "formula") && length(random) == 2L,
        "random must be a one-sided formula, for example ~1 or ~1 + time.")
  .need(sum(!vapply(list(relationship, relationship_inverse), is.null, TRUE)) == 1L,
        "Supply exactly one of relationship or relationship_inverse.")
  .need(!is.null(genetic_variance) && !is.null(residual_variance),
        "Supply genetic_variance and residual_variance from the model; these are not estimated automatically.")
  animals <- as.character(data[[animal]])
  focal_ids <- as.character(focal_ids);.ids(focal_ids, length(focal_ids), "focal animal IDs")
  .need(length(focal_ids) > 0, "Select at least one focal animal.")
  {
    K <- if (!is.null(relationship_inverse)) relationship_inverse else relationship
    .matrix_input(K, "relationship matrix")
    .ids(rownames(K), nrow(K), "relationship animal IDs")
    K <- .named_square(K, rownames(K), "relationship matrix")
    K <- .sparse(.symmetric(.sparse(K), "relationship matrix", tolerance))
    source <- if (is.null(relationship_inverse)) "supplied relationship" else "supplied relationship inverse"
    if (is.null(relationship_inverse)) {
      .memory_guard(24 * as.double(nrow(K))^2, budget, "Relationship inversion")
      factor <- .spd_factor(K, "relationship")
      inverse <- Matrix::solve(factor, Matrix::Diagonal(nrow(K)), system = "A")
      dimnames(inverse) <- dimnames(K);K <- inverse
    }
  }
  ids <- rownames(K);N <- length(ids)
  frame <- stats::model.frame(formula, data, na.action = stats::na.pass, drop.unused.levels = TRUE)
  .need(is.null(stats::model.offset(frame)), "Formula offsets are not supported in this version.")
  response <- stats::model.response(frame)
  .need(is.numeric(response) && !is.complex(response), "The response must be numeric.")
  if (is.null(dim(response))) response <- matrix(response, ncol = 1,
    dimnames = list(NULL, paste(deparse(formula[[2]]), collapse = "")))
  .need(nrow(response) == nrow(data) && !any(is.infinite(response)), "Invalid response dimensions or infinite values.")
  nt <- ncol(response);trait_names <- colnames(response)
  if (is.null(trait_names)) trait_names <- paste0("trait", seq_len(nt))
  .ids(trait_names, nt, "response names")
  observed <- !is.na(response);.need(all(colSums(observed) > 0), "Every response needs observed records.")
  relevant <- rowSums(observed) > 0
  missing_id <- is.na(animals) | !nzchar(animals)
  .need(allow_missing || !any(missing_id[relevant]),
        "Observed random-effect IDs contain missing or empty values; explicitly use allow_missing=TRUE only for zero incidence.")
  .need(all(c(animals[relevant & !missing_id], focal_ids) %in% ids),
        "Some observed or focal animals are absent from the relationship matrix.")
  basis_frame <- stats::model.frame(random, data[relevant, , drop = FALSE], na.action = stats::na.fail)
  basis_terms <- stats::terms(basis_frame)
  basis <- stats::model.matrix(basis_terms, basis_frame)
  basis_contrasts <- attr(basis, "contrasts")
  .need(ncol(basis) > 0 && all(is.finite(basis)), "Random-effect basis must contain finite columns.")
  .need(nt == 1L || (ncol(basis) == 1L && all(basis == 1)),
        "Multiple responses currently support an animal intercept (~1); multi-trait random regression is not yet supported.")
  k <- if (nt > 1) nt else ncol(basis)
  components <- if (nt > 1) trait_names else colnames(basis)
  G <- genetic_variance
  if (is.numeric(G) && is.null(dim(G)) && length(G) == 1L && k == 1L) G <- matrix(G, 1, 1)
  .matrix_input(G, "genetic_variance")
  .need(identical(dim(G), c(k, k)), "genetic_variance dimensions must match the random components.")
  if (!is.null(rownames(G)) || !is.null(colnames(G))) G <- .named_square(G, components, "genetic_variance")
  G <- as.matrix(.symmetric(G, "genetic_variance", tolerance))
  .need(isTRUE(tryCatch({chol(G); TRUE}, error = function(e) FALSE)), "genetic_variance must be positive definite.")
  K <- .sparse(K)
  .memory_guard(24 * as.double(length(K@x)) * k^2 + 32 * (as.double(length(focal_ids)) * k)^2,
                budget, "Joint random precision preparation")
  .need(as.double(N) * k < .Machine$integer.max && as.double(length(K@x)) * k^2 < .Machine$integer.max,
        "Joint precision exceeds the current index bound.")
  coef_ids <- if (k == 1) ids else make.unique(paste(rep(components, each = N), rep(ids, k), sep = ":"))
  Q <- NULL
  if (prepare_precision) {
    Q <- Matrix::kronecker(Matrix::Matrix(solve(G), sparse = TRUE), K)
    dimnames(Q) <- list(coef_ids, coef_ids)
  }
  fixed_terms <- stats::delete.response(stats::terms(frame));pieces <- list();removed <- list()
  if (is.null(template)) for (t in seq_len(nt)) {
    rows <- which(observed[, t])
    fixed_frame <- stats::model.frame(fixed_terms, data[rows, , drop = FALSE],
                                      na.action = stats::na.fail, drop.unused.levels = FALSE)
    for (j in seq_along(fixed_frame)) {
      if (is.character(fixed_frame[[j]])) fixed_frame[[j]] <- factor(fixed_frame[[j]])
      if (is.factor(fixed_frame[[j]]) && nlevels(fixed_frame[[j]]) < 2L)
        fixed_frame[[j]] <- rep(1, nrow(fixed_frame))
    }
    Xpart <- Matrix::sparse.model.matrix(fixed_terms, fixed_frame)
    independent <- .identified_design(Xpart, tolerance)
    pieces[[t]] <- independent$X;removed[[t]] <- independent$removed
  }
  locations <- which(observed, arr.ind = TRUE) # trait-major stacking
  nr <- nrow(locations);row <- locations[, 1];trait <- locations[, 2]
  X <- if (!is.null(template)) template$inputs$X else
    if (nt == 1) pieces[[1]] else do.call(Matrix::bdiag, pieces)
  if (!is.null(template)) {
    X <- template$inputs$X;removed <- template$summary$redundant_fixed_columns
  }
  known <- which(!missing_id[row])
  if (nt == 1) {
    basis <- basis[match(row, which(relevant)), , drop = FALSE]
    Z <- Matrix::sparseMatrix(i = rep(known, k),
      j = rep(match(animals[row[known]], ids), k) + rep((seq_len(k) - 1L) * N, each = length(known)),
      x = as.numeric(basis[known, , drop = FALSE]), dims = c(nr, N * k))
  } else Z <- Matrix::sparseMatrix(i = known,
      j = (trait[known] - 1L) * N + match(animals[row[known]], ids), x = 1, dims = c(nr, N * k))
  record_ids <- paste0("row", row, ":trait", trait)
  rownames(X) <- rownames(Z) <- record_ids;colnames(Z) <- coef_ids
  if (!is.null(template)) W <- template$inputs$residual_precision else {
  rv <- residual_variance
  if (inherits(rv, "formula")) {
    .need(nt == 1L && length(rv) == 2L, "A residual variance formula is supported for a single response only.")
    rv <- eval(rv[[2]], data, environment(rv))
  }
  if (nt == 1L) {
    .need(is.numeric(rv) && is.null(dim(rv)) && length(rv) %in% c(1L, nrow(data)),
          "residual_variance must be a positive scalar or one variance per original row.")
    rv <- rep_len(rv, nrow(data))[row]
    .need(all(is.finite(rv)) && all(rv > 0), "Observed residual variances must be finite and positive.")
    W <- 1 / rv
  } else {
    if (is.numeric(rv) && is.null(dim(rv))) {
      .need(length(rv) %in% c(1L, nt), "Supply a residual variance for each response or a residual covariance matrix.")
      if (!is.null(names(rv))) {
        .ids(names(rv), nt, "residual variance names")
        .need(setequal(names(rv), trait_names), "Residual variance trait names differ.")
        rv <- rv[trait_names]
      }
      rv <- diag(rep_len(rv, nt))
    }
    .matrix_input(rv, "residual_variance")
    .need(identical(dim(rv), c(nt, nt)), "Residual covariance dimensions do not match the responses.")
    if (!is.null(rownames(rv)) || !is.null(colnames(rv))) rv <- .named_square(rv, trait_names, "residual_variance")
    rv <- as.matrix(.symmetric(rv, "residual_variance", tolerance))
    .need(isTRUE(tryCatch({chol(rv); TRUE}, error = function(e) FALSE)), "Residual covariance must be positive definite.")
    if (all(rv[row(rv) != col(rv)] == 0) || all(rowSums(observed) <= 1)) W <- 1 / diag(rv)[trait] else {
      chunks <- split(seq_len(nr), row);ii <- jj <- numeric();xx <- numeric()
      for (positions in chunks) {
        block <- solve(rv[trait[positions], trait[positions], drop = FALSE])
        ii <- c(ii, rep(positions, times = length(positions)))
        jj <- c(jj, rep(positions, each = length(positions)));xx <- c(xx, as.numeric(block))
      }
      W <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(nr, nr),
                               dimnames = list(record_ids, record_ids))
    }
  }
  }
  focal_index <- unlist(lapply(seq_len(k), function(c) (c - 1L) * N + match(focal_ids, ids)))
  kind <- if (nt > 1) "multiple responses" else if (k > 1 || any(basis != 1)) "random regression" else "animal model"
  list(term = list(inverse = K, covariance = G, ids = ids,
                   components = components, zero_incidence_rows = sum(missing_id[relevant])),
       inputs = list(X = X, Z = Z, random_precision = Q, focal_ids = coef_ids[focal_index],
                     residual_precision = W, y = response[locations]),
       targets = list(kind = kind, animal_ids = focal_ids, components = components,
         basis_terms = basis_terms, basis_contrasts = basis_contrasts,
         basis_xlevels = lapply(basis_frame[vapply(basis_frame, is.factor, TRUE)], levels)),
       summary = list(formula = formula, random = random, relationship = source,
         animals = N, focal_animals = length(focal_ids), observed_values = nr,
         excluded_response_values = sum(!observed), components = components,
         redundant_fixed_columns = removed))
}
