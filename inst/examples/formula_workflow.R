# All relationships, IDs and observations are synthetic.
library(varestcf)
source(system.file("examples", "synthetic_models.R", package = "varestcf"))
models <- synthetic_models()
# The fixture supplies a previously constructed inverse relationship matrix.
Ainv <- models$animal_pedigree_relationship$random_precision * 4
records <- data.frame(
  animal = rep(rownames(Ainv), each = 3),
  group = rep(c(-.5, .5), 12),
  time = seq(-1.2, 1.2, length.out = 24),
  weight = as.numeric(models$animal_pedigree_relationship$y)
)
focals <- rownames(Ainv)[c(3, 5, 6, 8)]
fit <- reconstruct_focal_pev(
  data = records, formula = weight ~ group, animal = "animal",
  relationship_inverse = Ainv, focal_ids = focals,
  genetic_variance = 4, residual_variance = 5
)
print(fit)
print(varest_cf(fit))

Hinv <- models$animal_combined_relationship$random_precision * 4
fit_rr <- reconstruct_focal_pev(
  data = records, formula = weight ~ group, animal = "animal",
  relationship_inverse = Hinv, focal_ids = focals,
  random = ~1 + time,
  genetic_variance = matrix(c(4, .3, .3, .7), 2),
  residual_variance = ~exp(log(5) + .25 * time)
)
print(varest_cf(fit_rr, at = c(-1, 0, 1)))