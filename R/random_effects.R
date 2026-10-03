# Multiple terms share fixed designs/residuals; only their random designs differ.
# Independent terms form block-diagonal precision. Correlated terms use the
# inverse of the JOINT component covariance, never marginal inverses plus a patch.
.component_covariance <- function(x, labels, name, tolerance) {
  k <- length(labels)
  if (is.numeric(x) && is.null(dim(x)) && length(x) == 1L && k == 1L) x <- matrix(x, 1, 1)
  .matrix_input(x, name)
  .need(identical(dim(x), c(k, k)), paste0(name, " dimensions must match its components."))
  if (!is.null(rownames(x)) || !is.null(colnames(x))) x <- .named_square(x, labels, name)
  x <- as.matrix(.symmetric(x, name, tolerance))
  .need(isTRUE(tryCatch({chol(x); TRUE}, error = function(e) FALSE)),
        paste0(name, " must be positive definite."))
  x
}

.formula_effects <- function(data, formula, focal_ids, random_effects, focal_effect,
                             correlated_effects, residual_variance, budget, tolerance) {
  .need(is.list(random_effects) && length(random_effects) > 0L,
        "random_effects must be a nonempty named list.")
  effect_names <- names(random_effects)
  .ids(effect_names, length(random_effects), "random_effects names")
  .need(is.character(focal_effect) && length(focal_effect) == 1L &&
          !is.na(focal_effect) && focal_effect %in% effect_names,
        "Select one focal_effect by its name in random_effects.")
  .need(is.data.frame(data) && nrow(data) > 0L, "data must be a nonempty data.frame.")
  .need(inherits(formula, "formula") && length(formula) == 3L, "formula must include a response.")
  .need(!is.null(residual_variance), "Supply residual_variance from the model.")
  response <- stats::model.response(stats::model.frame(formula, data, na.action = stats::na.pass))
  .need(is.numeric(response) && !is.complex(response), "The response must be numeric.")
  relevant <- if (is.null(dim(response))) !is.na(response) else rowSums(!is.na(response)) > 0
  .need(any(relevant), "The model needs observed responses.")
  nt <- if (is.null(dim(response))) 1L else ncol(response)
  prepared <- vector("list", length(effect_names));names(prepared) <- effect_names
  supplied_variance <- setNames(logical(length(effect_names)), effect_names)
  template <- NULL
  # Prepare the focal term first, so other terms reuse its fixed/residual inputs.
  for (name in c(focal_effect, setdiff(effect_names, focal_effect))) {
    spec <- random_effects[[name]]
    .need(is.list(spec) && !is.null(names(spec)) && !anyDuplicated(names(spec)) &&
            all(names(spec) %in% c("id", "random", "variance", "relationship",
              "relationship_inverse", "levels", "allow_missing")),
          paste0("Invalid fields in random effect ", name, "."))
    .need(is.character(spec$id) && length(spec$id) == 1L && spec$id %in% names(data),
          paste0("Effect ", name, " needs id naming a data column."))
    allow_missing <- if (is.null(spec$allow_missing)) FALSE else spec$allow_missing
    .need(is.logical(allow_missing) && length(allow_missing) == 1L && !is.na(allow_missing),
          "allow_missing must be TRUE or FALSE.")
    random <- if (is.null(spec$random)) ~1 else spec$random
    .need(inherits(random, "formula") && length(random) == 2L,
          paste0("Effect ", name, " needs a one-sided random formula."))
    relation_count <- sum(!vapply(list(spec[["relationship"]], spec[["relationship_inverse"]]), is.null, TRUE))
    .need(relation_count <= 1L, paste0("Effect ", name, ": supply exactly one relationship input."))
    .need(name != focal_effect || relation_count == 1L,
          "The focal genetic effect requires a supplied relationship or relationship_inverse (identity is allowed explicitly).")
    independent <- relation_count == 0L
    if (independent) {
      observed_ids <- as.character(data[[spec$id]][relevant])
      observed_ids <- observed_ids[!is.na(observed_ids) & nzchar(observed_ids)]
      levels <- if (is.null(spec$levels)) unique(observed_ids) else as.character(spec$levels)
      .ids(levels, length(levels), paste0(name, " levels"))
      .need(length(levels) > 0L, paste0("Effect ", name, " needs at least one level."))
      inverse <- Matrix::Diagonal(length(levels))
      dimnames(inverse) <- list(levels, levels)
    } else {
      .need(is.null(spec$levels), "levels is only used for an independent effect without a relationship matrix.")
      inverse <- spec[["relationship_inverse"]]
      relation <- if (is.null(inverse)) spec[["relationship"]] else inverse
      .matrix_input(relation, paste0(name, " relationship"))
      levels <- rownames(relation);.ids(levels, nrow(relation), paste0(name, " relationship IDs"))
    }
    variance <- spec$variance;supplied_variance[name] <- !is.null(variance)
    if (is.null(variance)) {
      k <- if (nt > 1L) nt else ncol(stats::model.matrix(random,
        stats::model.frame(random, data[relevant, , drop = FALSE], na.action = stats::na.fail)))
      variance <- diag(k) # temporary marginal; the joint group covariance is mandatory below
    }
    fit <- .formula_model(data, formula, spec$id,
      if (name == focal_effect) focal_ids else levels[1L],
      if (!independent) spec[["relationship"]] else NULL, inverse, random,
      variance, residual_variance, budget, tolerance,
      allow_missing = allow_missing, template = template, prepare_precision = FALSE)
    if (independent) fit$summary$relationship <- "independent levels"
    prepared[[name]] <- fit
    if (is.null(template)) template <- fit
  }
  counts <- vapply(prepared, function(z) ncol(z$inputs$Z), 1L)
  offsets <- c(0, utils::head(cumsum(as.double(counts)), -1L))
  names(offsets) <- effect_names
  .need(sum(as.double(counts)) < .Machine$integer.max, "Joint precision exceeds the current index bound.")
  groups <- if (is.null(correlated_effects)) list() else correlated_effects
  .need(is.list(groups), "correlated_effects must be a list of joint covariance groups.")
  used <- character();units <- list()
  for (group in groups) {
    .need(is.list(group) && setequal(names(group), c("effects", "covariance")) &&
            !anyDuplicated(names(group)), "Each correlated group needs effects and covariance only.")
    names_group <- group$effects
    .need(is.character(names_group) && length(names_group) >= 2L && !anyNA(names_group) &&
            !anyDuplicated(names_group) && all(names_group %in% effect_names) &&
            !any(names_group %in% used), "Correlated groups need distinct existing effects; groups cannot overlap.")
    labels <- unlist(lapply(names_group, function(name) {
      components <- prepared[[name]]$term$components
      if (length(components) == 1L) name else paste(name, components, sep = ":")
    }), use.names = FALSE)
    .ids(labels, length(labels), "Correlated component labels")
    covariance <- .component_covariance(group$covariance, labels, "Joint covariance", tolerance)
    reference <- prepared[[names_group[1L]]]$term
    position <- 0L
    for (name in names_group) {
      term <- prepared[[name]]$term
      .need(setequal(term$ids, reference$ids), "Correlated effects must share the same relationship IDs.")
      aligned <- term$inverse[reference$ids, reference$ids, drop = FALSE]
      .need(max(abs(aligned-reference$inverse)) <= tolerance *
              max(abs(reference$inverse), .Machine$double.xmin),
            "Correlated effects must share the same relationship inverse; use matrix inputs for a general cross-kernel.")
      k <- length(term$components);idx <- position + seq_len(k)
      if (supplied_variance[name])
        .need(max(abs(term$covariance-covariance[idx, idx, drop = FALSE])) <= tolerance *
                max(abs(term$covariance), .Machine$double.xmin),
              "A term variance conflicts with its joint covariance diagonal block.")
      position <- position + k
    }
    # Unit columns share reference-ID order; map them back to each term's own order.
    positions <- unlist(lapply(names_group, function(name) {
      term <- prepared[[name]]$term
      offsets[name] + unlist(lapply(seq_along(term$components), function(c)
        (c-1L)*length(term$ids)+match(reference$ids, term$ids)), use.names = FALSE)
    }), use.names = FALSE)
    units[[length(units)+1L]] <- list(inverse=reference$inverse, covariance=covariance, positions=positions)
    used <- c(used, names_group)
  }
  for (name in setdiff(effect_names, used)) {
    .need(supplied_variance[name], paste0("Supply variance for effect ", name, " or include it in a correlated group."))
    term <- prepared[[name]]$term
    units[[length(units)+1L]] <- list(inverse=term$inverse, covariance=term$covariance,
                                    positions=offsets[name]+seq_len(counts[name]))
  }
  bound <- sum(vapply(units, function(unit)
    as.double(length(unit$inverse@x))*nrow(unit$covariance)^2, 0))
  .need(bound < .Machine$integer.max, "Joint precision exceeds the current index bound.")
  input_bytes <- sum(vapply(prepared, function(z)
    as.numeric(object.size(z$inputs$Z))+as.numeric(object.size(z$term$inverse)), 0))
  .memory_guard(input_bytes + as.numeric(object.size(template$inputs$X)) +
      as.numeric(object.size(template$inputs$residual_precision)) +
      32*as.double(length(template$inputs$focal_ids))^2+48*bound, budget,
      "Multiple random terms and joint precision preparation")
  blocks <- lapply(units, function(unit)
    Matrix::kronecker(Matrix::Matrix(solve(unit$covariance), sparse=TRUE), unit$inverse))
  unit_positions <- unlist(lapply(units, function(unit) unit$positions), use.names=FALSE)
  Q <- do.call(Matrix::bdiag, blocks)
  permutation <- order(unit_positions);Q <- Q[permutation, permutation, drop=FALSE]
  Z <- do.call(cbind, lapply(prepared, function(z) z$inputs$Z))
  coefficient_ids <- make.unique(unlist(lapply(effect_names, function(name)
    paste(name, colnames(prepared[[name]]$inputs$Z), sep=":")), use.names=FALSE))
  dimnames(Q) <- list(coefficient_ids, coefficient_ids);colnames(Z) <- coefficient_ids
  focal <- prepared[[focal_effect]]
  focal_index <- offsets[focal_effect] +
    match(focal$inputs$focal_ids, colnames(focal$inputs$Z))
  focal$inputs$Z <- Z;focal$inputs$random_precision <- Q
  focal$inputs$focal_ids <- coefficient_ids[focal_index]
  focal$summary$focal_effect <- focal_effect
  focal$summary$random_effects <- lapply(prepared, function(z) list(
    levels=length(z$term$ids), components=z$term$components,
    relationship=z$summary$relationship, zero_incidence_rows=z$term$zero_incidence_rows))
  focal$summary$correlated_effects <- lapply(groups, function(group) group$effects)
  focal$term <- NULL
  focal
}
