# Entirely invented animals, pedigree, relationships and observations.
# No genotypes or identifiers from the private data are used.
synthetic_models <- function() {
  ids <- paste0("animal", seq_len(8))
  sire <- c(0L, 0L, 1L, 1L, 3L, 3L, 5L, 5L)
  dam  <- c(0L, 0L, 2L, 2L, 4L, 4L, 6L, 6L)
  A <- matrix(0, 8, 8, dimnames = list(ids, ids))
  for (i in seq_len(8)) {
    if (i > 1L) for (j in seq_len(i - 1L)) {
      A[i, j] <- (if (sire[i]) A[sire[i], j] else 0) / 2 +
                 (if (dam[i]) A[dam[i], j] else 0) / 2
      A[j, i] <- A[i, j]
    }
    A[i, i] <- 1 + if (sire[i] && dam[i]) A[sire[i], dam[i]] / 2 else 0
  }
  Ai <- solve(A)
  genotyped <- c(3L, 5L, 6L, 8L)
  feature <- matrix(sin(seq_len(20)), 4, 5)
  Grel <- .8 * A[genotyped, genotyped] + .2 * tcrossprod(feature) / 5 + diag(.1, 4)
  Hi <- Ai
  Hi[genotyped, genotyped] <- Hi[genotyped, genotyped] + solve(Grel) - solve(A[genotyped, genotyped])
  animal <- rep(seq_len(8), each = 3)
  nr <- length(animal); record <- paste0("record", seq_len(nr))
  x <- seq(-1.2, 1.2, length.out = nr)
  X <- cbind(intercept = 1, group = rep(c(-.5, .5), nr / 2))
  rownames(X) <- record
  incidence <- diag(8)[animal, , drop = FALSE]
  dimnames(incidence) <- list(record, ids)
  y <- setNames(30 + X[, 2] + sin(seq_len(nr) * 1.7), record)
  make <- function(XX, ZZ, Q, w, focal) {
    dimnames(Q) <- list(colnames(ZZ), colnames(ZZ))
    list(X = XX, Z = ZZ, random_precision = Q, residual_precision = w,
         focal_ids = focal, y = y)
  }
  out <- list(animal_pedigree_relationship = make(X, incidence, Ai / 4, .2, ids[c(8, 1, 5, 3)]),
              animal_combined_relationship = make(X, incidence, Hi / 4, .2, ids[c(8, 1, 5, 3)]))
  coef_ids <- c(paste0("a0:", ids), paste0("a1:", ids))
  Zrr <- cbind(incidence, incidence * x);colnames(Zrr) <- coef_ids
  G0 <- matrix(c(4, .3, .3, .7), 2)
  Qrr <- kronecker(solve(G0), Hi)
  focal <- c(paste0("a0:", ids[genotyped]), paste0("a1:", ids[genotyped]))
  out$random_regression_homogeneous <- make(X, Zrr, Qrr, .2, focal)
  out$random_regression_heterogeneous <- make(X, Zrr, Qrr, exp(-log(5) - .25 * x), focal)
  trait <- rep(seq_len(4), length.out = nr)
  Zmulti <- matrix(0, nr, 32)
  Zmulti[cbind(seq_len(nr), (trait - 1L) * 8L + animal)] <- 1
  dimnames(Zmulti) <- list(record, paste0(rep(paste0("trait", 1:4), each = 8), ":", rep(ids, 4)))
  Xmulti <- diag(4)[trait, , drop = FALSE];rownames(Xmulti) <- record
  Gmulti <- diag(c(4, 5, 3, 6)) + .3
  out$multivariate_animal <- make(Xmulti, Zmulti, kronecker(solve(Gmulti), Hi), 1 / c(5, 6, 7, 8)[trait],
                 colnames(Zmulti)[unlist(lapply(0:3, function(k) 8L * k + genotyped))])
  # Multiple random terms: focal genetic effects remain a subset of the full model.
  repeated_Z <- cbind(incidence, incidence)
  colnames(repeated_Z) <- c(paste0("genetic:", ids), paste0("permanent:", ids))
  repeated_Q <- as.matrix(Matrix::bdiag(Ai / 4, diag(8) / 2))
  out$repeatability <- make(X, repeated_Z, repeated_Q, .2, paste0("genetic:", ids[c(8, 1, 5, 3)]))
  maternal_Z <- matrix(0, nr, 8)
  observed_dam <- dam[animal];known_dam <- which(observed_dam > 0)
  maternal_Z[cbind(known_dam, observed_dam[known_dam])] <- 1
  joint_Z <- cbind(incidence, maternal_Z)
  colnames(joint_Z) <- c(paste0("direct:", ids), paste0("maternal:", ids))
  direct_maternal_G <- matrix(c(4, .4, .4, 2), 2)
  out$direct_maternal <- make(X, joint_Z, kronecker(solve(direct_maternal_G), Ai),
    .2, paste0("direct:", ids[c(8, 1, 5, 3)]))
  out
}

# Only for tiny synthetic examples and independent test oracles: dense full inverse.
dense_reference <- function(model) {
  D <- cbind(model$X, model$Z);p <- ncol(model$X);r <- ncol(model$Z)
  W <- model$residual_precision
  if (is.null(dim(W))) W <- diag(rep_len(W, nrow(D)))
  C <- crossprod(D, W %*% D)
  C[p + seq_len(r), p + seq_len(r)] <- C[p + seq_len(r), p + seq_len(r)] + model$random_precision
  inverse <- solve(C)
  u <- as.numeric(solve(C, crossprod(D, W %*% model$y)))[p + seq_len(r)]
  names(u) <- colnames(model$Z)
  idx <- p + match(model$focal_ids, colnames(model$Z))
  S <- inverse[idx, idx, drop = FALSE];dimnames(S) <- list(model$focal_ids, model$focal_ids)
  list(Sigma = S, ebv = u[model$focal_ids], C = C, p = p, inverse = inverse)
}
