# Entirely synthetic user workflow for multiple random effects.
library(varestcf)
source(system.file("examples", "synthetic_models.R", package = "varestcf"))
models <- synthetic_models()
Kinv <- models$animal_pedigree_relationship$random_precision*4
records <- data.frame(
  animal=rep(rownames(Kinv), each=3),
  group=rep(c(-.5,.5),12),
  trait=as.numeric(models$animal_pedigree_relationship$y),
  dam=rep(c(NA,NA,"animal2","animal2","animal4","animal4","animal6","animal6"),each=3)
)
focal_ids <- rownames(Kinv)[c(8,1,5,3)]

fit_repeatability <- reconstruct_focal_pev(
  data=records, formula=trait~group,
  random_effects=list(
    additive=list(id="animal",relationship_inverse=Kinv,variance=4),
    permanent=list(id="animal",variance=2)
  ),
  focal_effect="additive",focal_ids=focal_ids,residual_variance=5
)
print(fit_repeatability)
print(varest_cf(fit_repeatability))

# In THIS synthetic model, founder records have no maternal term.
# Unknown parental effects in another model must instead be represented explicitly.
Gdm <- matrix(c(4,.4,.4,2),2,
              dimnames=list(c("direct","maternal"),c("direct","maternal")))
fit_maternal <- reconstruct_focal_pev(
  data=records,formula=trait~group,
  random_effects=list(
    direct=list(id="animal",relationship_inverse=Kinv),
    maternal=list(id="dam",relationship_inverse=Kinv,allow_missing=TRUE)
  ),
  correlated_effects=list(list(effects=c("direct","maternal"),covariance=Gdm)),
  focal_effect="direct",focal_ids=focal_ids,residual_variance=5
)
print(fit_maternal)
print(varest_cf(fit_maternal))
