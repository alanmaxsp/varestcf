# Entirely synthetic data and a supplied relationship matrix.
library(varestcf)

ids <- paste0("animal", 1:4)
K <- matrix(c(
  1,   0,   0.5, 0.5,
  0,   1,   0.5, 0.5,
  0.5, 0.5, 1,   0.5,
  0.5, 0.5, 0.5, 1
), nrow = 4, byrow = TRUE, dimnames = list(ids, ids))

records <- data.frame(
  animal = rep(ids[1:3], each = 2),
  group = rep(c("group1", "group2"), 3),
  trait = c(10, 11, 9, 10, 13, 12)
)

fit <- reconstruct_focal_pev(
  data = records,
  formula = trait ~ group,
  animal = "animal",
  relationship = K,
  focal_ids = ids[2:4],
  genetic_variance = 2,
  residual_variance = 1
)

result <- varest_cf(fit)
print(result)

# Compatible focal breeding values and their complete error covariance.
print(fit$ebv)
print(fit$Sigma)
