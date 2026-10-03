\name{reconstruct_focal_pev}
\alias{reconstruct_focal_pev}
\title{Reconstruct a complete focal prediction error covariance block}
\usage{
reconstruct_focal_pev(X = NULL, Z = NULL, random_precision = NULL, focal_ids,
  residual_precision = 1, block_size = 8L,
  memory_budget_gib = 1, tolerance = 1e-8,
  data = NULL, formula = NULL, animal = NULL,
  relationship = NULL, relationship_inverse = NULL, random = ~1,
  genetic_variance = NULL, residual_variance = NULL, y = NULL,
  random_effects = NULL, focal_effect = NULL, correlated_effects = NULL)
}
\arguments{
\item{X}{Numeric matrix or double Matrix object of identified fixed-effect
  covariates. Columns must be linearly independent. Zero columns are allowed.}
\item{Z}{Numeric matrix or double Matrix object of random-effect covariates.
  Column names are unique coefficient identifiers. Row order must match X.}
\item{random_precision}{Named, symmetric positive definite matrix Q equal to
  the inverse covariance of all random coefficients in Z, including variance
  scale. Row and column names are aligned to Z column names. This is not merely
  a relationship inverse unless its variance scale is one.}
\item{focal_ids}{Unique animal IDs in formula mode, or random coefficients in matrix mode,
  in the desired output order. Coefficients can represent multiple traits
  or random regression terms.}
\item{residual_precision}{Inverse residual covariance. A positive scalar,
  a positive vector with one entry per record, or a named symmetric positive
  definite matrix. Named vectors/matrices require matching row names in X/Z;
  their entries are aligned by record identifier. Unnamed vectors are positional.}
\item{block_size}{Positive integer: columns of the inverse solved together.}
\item{memory_budget_gib}{Positive advisory working-storage budget in GiB.
  It is not an operating-system memory limit or a guarantee of sufficient RAM.}
\item{tolerance}{Positive relative tolerance for symmetry and normalized
  backward errors of the selected-column solves.}
\item{data}{Data frame for the automatic formula interface.}
\item{formula}{Response and fixed effects, e.g. weight ~ sex + group.
  Use cbind(y1,y2) for multiple numeric responses. Missing responses are allowed.}
\item{animal}{Name of the animal-ID column in data.}
\item{relationship}{Named positive definite relationship matrix A, H or another
  valid relationship for all animals in the model, without genetic variance scaling.
  Supply exactly one of relationship or relationship_inverse in formula mode.}
\item{relationship_inverse}{Precomputed inverse of the relationship matrix.
  Preferred to avoid an additional, potentially dense inversion. All animal IDs
  must appear in both row and column names. No pedigree or genotypes are processed.}
\item{random}{One-sided formula for the animal's random coefficient basis;
  default ~1. Random regression, e.g. ~1 + time, supports a single response.}
\item{genetic_variance}{Known scalar genetic variance or covariance matrix
  between random coefficients/responses. Named matrices are aligned to components.}
\item{residual_variance}{For one response: known positive scalar, vector per original
  data row, or one-sided expression evaluated in data. For multiple responses:
  scalar, vector per response, or residual covariance matrix.}
\item{y}{Optional finite response vector in explicit matrix mode, in design row
  order. Enables calculation of focal EBV from the same equations.}
\item{random_effects}{Alternative named list of random term specifications.
  Each term requires \code{id} (a data column). Optional fields are
  \code{random} (one-sided basis, default ~1), \code{variance} (scalar or
  component covariance), \code{relationship} or \code{relationship_inverse},
  \code{levels}, and \code{allow_missing}. Without a relationship input,
  levels are independent and an identity relationship is built; optional levels
  give their complete order, otherwise observed levels are inferred. The focal
  genetic term requires an explicitly supplied relationship input.
  Do not mix this specification with animal, random, genetic_variance or the
  top-level relationship inputs.}
\item{focal_effect}{Required name of the genetic term to summarize when
  random_effects is supplied. All other terms remain in the complete equations.}
\item{correlated_effects}{Optional list of disjoint groups, each with
  \code{effects} (at least two term names) and \code{covariance} (joint component
  covariance). Group terms must have the same relationship inverse and ID set.
  The covariance follows the effects order, then each term's component order.
  A scalar term label is its name; multi-component labels are term:component.
  Named covariance matrices are aligned by these labels. Within a group, a
  term variance may be omitted; if supplied, it must match the joint diagonal
  block. Outside groups every variance is required.}
}
\description{
Forms the complete identified mixed-model equations, factors them once, and
solves only the columns needed for the focal covariance block. Fixed-effect
uncertainty and cross-coefficient covariances are retained.
The data/formula interface constructs designs and precision matrices and obtains
compatible EBV automatically from the supplied relationships and variance components.
}
\details{
For a model with fixed design X, random design Z, residual precision W and
random precision Q, the coefficient matrix has blocks X'WX, X'WZ, Z'WX,
and Z'WZ + Q. If E selects the focal random coefficients in this full system,
the result is E'C-inverse E. The function does not invert the focal principal
submatrix of C, and does not treat fixed effects as known.

The covariance is conditional on the supplied model and variance components.
It does not include uncertainty of estimated components. An equivalent
identified parameterization of X must be supplied in matrix mode; redundant columns
are rejected. Formula mode removes and reports redundant fixed columns using a
scaled sparse QR decomposition. No variance components are estimated.

For coefficient order grouped by component, Q can be constructed as
kronecker(solve(G0), K-inverse). K can be A, H or another valid covariance.
Formula mode constructs Q from the supplied relationship inverse and genetic
covariance; matrix mode receives Q explicitly. A is never substituted for H.
Relationships are supplied for all model animals, not just focal animals.

Multiple responses have separate fixed-effect coefficients and trait-major stacking.
With missing responses and correlated residuals, each observed residual covariance
submatrix is inverted separately. Multiple random terms can be independent or
correlated. Independent terms form separate precision blocks; for a correlated
group the joint component covariance is inverted before taking its Kronecker
product with the common relationship inverse. Additional terms may have their
own random-regression bases for a single response.

Missing or empty term IDs are rejected on observed rows by default.
Setting allow_missing=TRUE explicitly assigns zero incidence for that term;
it does not infer or model an unknown parent's effect. Represent an unknown
level explicitly when that effect is present in the intended model.
Unobserved levels in supplied relationships or explicit levels are retained.

Formula offsets, multi-trait random regression, general cross-kernel
covariances and heterogeneous multivariate residual blocks are not automated.
The scope is linear Gaussian mixed models conditional on supplied covariances.

The reference performs positive-definiteness checks of Q and of a general
residual precision by additional factorizations. Conservative assembly bounds
and 32-bit CSC index limits can reject a model even when a more optimized
representation would fit. Factor fill and transient allocations remain
uncertain; a budget check after factorization cannot prevent all out-of-memory
conditions. Sigma is dense and its final Cholesky validation costs cubic time
in the number of focal coefficients. No cache or implicit approximation is used.
}
\value{
A \code{varestcf_fit} list with \code{Sigma}, the dense named focal covariance;
\code{ebv}, compatible breeding values when responses were supplied;
formula model metadata when applicable;
including the focal term, model terms and explicitly zero incidence counts;
\code{focal_ids}, its order; and \code{diagnostics}, including dimensions,
stored entries, object sizes, an assembly bound, per-block absolute residuals
and backward errors, and elapsed seconds per stage. Object sizes and estimates
are not measured peak resident memory. Backward errors use the maximum
absolute residual divided by the infinity norm of C times the maximum
absolute solution entry plus one. They do not bound forward error for an
ill-conditioned model.
}
\seealso{\code{\link{varest_cf}}}
\examples{
ids <- c("a", "b", "c")
X <- matrix(1, 3, 1)
Z <- diag(3); colnames(Z) <- ids
Q <- diag(3) / 2; dimnames(Q) <- list(ids, ids)
focal <- reconstruct_focal_pev(X, Z, Q, c("c", "a"))
focal$Sigma
}
