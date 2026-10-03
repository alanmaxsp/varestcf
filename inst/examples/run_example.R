library(varestcf)
source(system.file("examples", "synthetic_models.R", package = "varestcf"))
model <- synthetic_models()$animal_combined_relationship
inputs <- model[setdiff(names(model), "y")]
focal <- do.call(reconstruct_focal_pev, inputs)
# The dense solution is only a source of EBV for this tiny, fully synthetic example.
ebv <- dense_reference(model)$ebv
result <- varest_cf(ebv, focal$Sigma, normalization = "n")
print(result[c("estimate", "ebv_term", "error_term", "n", "normalization")])
print(focal$diagnostics[c("equations", "focal_coefficients", "C_bytes", "factor_bytes", "seconds")])
