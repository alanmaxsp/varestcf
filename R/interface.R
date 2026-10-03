reconstruct_focal_pev <- function(X = NULL, Z = NULL, random_precision = NULL, focal_ids,
    residual_precision = 1, block_size = 8L, memory_budget_gib = 1, tolerance = 1e-8,
    data = NULL, formula = NULL, animal = NULL,
    relationship = NULL, relationship_inverse = NULL, random = ~1,
    genetic_variance = NULL, residual_variance = NULL,
    y = NULL, random_effects = NULL, focal_effect = NULL, correlated_effects = NULL) {
  .scalar_positive(memory_budget_gib, "memory_budget_gib")
  .scalar_positive(tolerance, "tolerance")
  if (!is.null(data) || !is.null(formula)) {
    .need(is.null(X) && is.null(Z) && is.null(random_precision) && is.null(y) && missing(residual_precision),
          "Use either data/formula or explicit matrix inputs, not both.")
    prepared_at <- proc.time()[[3L]]
    if (is.null(random_effects)) {
      .need(is.null(focal_effect) && is.null(correlated_effects),
            "focal_effect/correlated_effects require random_effects.")
      prepared <- .formula_model(data, formula, animal, focal_ids,
        relationship, relationship_inverse, random, genetic_variance, residual_variance,
        memory_budget_gib, tolerance)
    } else {
      .need(is.null(animal) && is.null(relationship) && is.null(relationship_inverse) &&
              is.null(genetic_variance) && missing(random),
            "Put all random term inputs inside random_effects; do not mix with the single-effect shorthand.")
      prepared <- .formula_effects(data, formula, focal_ids, random_effects, focal_effect,
        correlated_effects, residual_variance, memory_budget_gib, tolerance)
    }
    preparation_seconds <- proc.time()[[3L]] - prepared_at
    out <- do.call(.reconstruct_from_matrices, c(prepared$inputs, list(
      block_size = block_size, memory_budget_gib = memory_budget_gib, tolerance = tolerance)))
    out$targets <- prepared$targets;out$model <- prepared$summary
    out$diagnostics$seconds <- c(preparation = preparation_seconds, out$diagnostics$seconds)
    out$diagnostics$seconds["total"] <- out$diagnostics$seconds["total"] + preparation_seconds
  } else {
    .need(is.null(relationship) && is.null(relationship_inverse) &&
            is.null(genetic_variance) && is.null(residual_variance) &&
            is.null(random_effects) && is.null(focal_effect) && is.null(correlated_effects),
          "Relationship and variance inputs require data and formula.")
    out <- .reconstruct_from_matrices(X, Z, random_precision, focal_ids,
      residual_precision, block_size, memory_budget_gib, tolerance, y)
  }
  class(out) <- "varestcf_fit"
  out
}

varest_cf <- function(ebv, Sigma = NULL, normalization = c("n", "n-1"), tolerance = 1e-8,
                      at = NULL, component = NULL) {
  normalization <- match.arg(normalization)
  if (!inherits(ebv, "varestcf_fit")) {
    .need(is.null(at) && is.null(component), "at/component are available for a reconstructed model object.")
    out <- .varest_from_covariance(ebv, Sigma, normalization, tolerance)
    class(out) <- "varestcf_result";return(out)
  }
  fit <- ebv
  .need(is.null(Sigma), "A model object already contains Sigma; do not supply a second matrix.")
  .need(!is.null(fit$ebv), "This matrix-only fit has no EBV. Supply y during reconstruction or pass EBV and Sigma explicitly.")
  .need(is.null(at) || is.null(component), "Choose at or component, not both.")
  if (is.null(fit$targets)) {
    .need(is.null(at) && is.null(component), "This matrix-only fit has no formula projection metadata.")
    return(varest_cf(fit$ebv, fit$Sigma, normalization, tolerance))
  }
  targets <- fit$targets;n <- length(targets$animal_ids);k <- length(targets$components)
  weights <- NULL;grid <- NULL
  if (!is.null(at)) {
    .need(targets$kind == "random regression", "at is used for random-regression models only.")
    variables <- all.vars(targets$basis_terms)
    if (is.numeric(at) && is.null(dim(at))) {
      .need(length(variables) == 1L, "Use at=data.frame(...) when the random formula has multiple variables.")
      at <- setNames(data.frame(at), variables)
    }
    .need(is.data.frame(at) && nrow(at) > 0, "at must be a nonempty data.frame or a numeric vector for one variable.")
    mf <- stats::model.frame(targets$basis_terms, at, na.action = stats::na.fail,
                            xlev = targets$basis_xlevels)
    weights <- stats::model.matrix(targets$basis_terms, mf, contrasts.arg = targets$basis_contrasts)
    .need(setequal(colnames(weights), targets$components) && all(is.finite(weights)),
          "Projection basis is incompatible with the fitted random effects.")
    weights <- weights[, targets$components, drop = FALSE];grid <- at
    labels <- apply(at, 1, function(z) paste(paste(names(at), z, sep = "="), collapse = ", "))
  } else {
    if (is.null(component)) {
      .need(targets$kind != "random regression",
            "For random regression use at=c(...) or at=data.frame(...); use component= to evaluate one coefficient.")
      component <- targets$components
    }
    .need(is.character(component) && length(component) > 0 && !anyDuplicated(component) &&
            all(component %in% targets$components), "Unknown or duplicated component.")
    weights <- diag(k)[match(component, targets$components), , drop = FALSE];labels <- component
  }
  results <- vector("list", nrow(weights))
  for (point in seq_len(nrow(weights))) {
    w <- weights[point, ];u <- setNames(numeric(n), targets$animal_ids)
    S <- matrix(0, n, n, dimnames = list(targets$animal_ids, targets$animal_ids))
    for (a in which(w != 0)) {
      ia <- (a - 1L) * n + seq_len(n);u <- u + w[a] * fit$ebv[ia]
      for (b in which(w != 0)) {
        ib <- (b - 1L) * n + seq_len(n)
        S <- S + w[a] * w[b] * fit$Sigma[ia, ib, drop = FALSE]
      }
    }
    names(u) <- targets$animal_ids;dimnames(S) <- list(targets$animal_ids, targets$animal_ids)
    result <- .varest_from_covariance(u, S, normalization, tolerance)
    results[[point]] <- data.frame(target = labels[point], n = n, estimate = result$estimate,
      ebv_term = result$ebv_term, error_term = result$error_term, normalization = normalization,
      row.names = NULL)
  }
  table <- do.call(rbind, results)
  structure(list(estimate = setNames(table$estimate, table$target), table = table,
    normalization = normalization, at = grid, diagnostics = list(source = "reconstructed model",
    focal_animals = targets$animal_ids, components = targets$components,
    focal_effect = fit$model$focal_effect)), class = "varestcf_result")
}

print.varestcf_fit <- function(x, ...) {
  cat("Complete focal prediction-error covariance\n")
  if (!is.null(x$model)) {
    cat("  Model:", x$targets$kind, "| Relationship:", x$model$relationship, "\n")
    cat("  Animals:", x$model$animals, "| Focals:", x$model$focal_animals,
        "| Observed values:", x$model$observed_values, "\n")
    cat("  Components:", paste(x$targets$components, collapse = ", "), "\n")
    if (!is.null(x$model$random_effects))
      cat("  Focal effect:", x$model$focal_effect, "| Model terms:",
          paste(names(x$model$random_effects), collapse = ", "), "\n")
  }
  cat("  Sigma:", nrow(x$Sigma), "x", ncol(x$Sigma), "| EBV:",
      if (is.null(x$ebv)) "not supplied" else "available from the same model", "\n")
  cat("  Elapsed:", round(x$diagnostics$seconds["total"], 3), "seconds\n")
  invisible(x)
}

print.varestcf_result <- function(x, ...) {
  cat("VarEst-CF\n")
  if (!is.null(x$diagnostics$focal_effect))
    cat("  Focal effect:", x$diagnostics$focal_effect, "\n")
  if (!is.null(x$table)) print(x$table, row.names = FALSE) else
    print(data.frame(n = x$n, estimate = x$estimate, ebv_term = x$ebv_term,
                     error_term = x$error_term, normalization = x$normalization), row.names = FALSE)
  invisible(x)
}
