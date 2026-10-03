\name{varest_cf}
\alias{varest_cf}
\title{Calculate focal genetic variance from EBV and complete covariance}
\usage{
varest_cf(ebv, Sigma = NULL, normalization = c("n", "n-1"), tolerance = 1e-8,
  at = NULL, component = NULL)
}
\arguments{
\item{ebv}{A reconstructed model object with EBV, or a finite named numeric
  vector of at least two focal breeding values.}
\item{Sigma}{Explicitly complete dense covariance matrix of prediction errors,
  with unique row and column identifiers matching the EBV names. Entries are
  aligned by identifier. Sparse Matrix inputs are rejected to avoid interpreting
  structural omissions as supplied zero covariances. Omit for a model object.}
\item{normalization}{Use denominator n (default) or n minus one.}
\item{tolerance}{Positive relative tolerance for covariance validation and
  numerical roundoff in the centered covariance trace.}
\item{at}{Random-regression evaluation points: numeric vector for one variable,
  or data frame whose columns match the random formula.}
\item{component}{Optional names of individual coefficients/responses to evaluate.
  Mutually exclusive with at. Random regression requires at or component.}
}
\description{
Calculates the sum of centered squared EBV and the centered prediction-error
covariance trace, divided by the chosen denominator. No mixed-model equations
or pedigree are reconstructed by this function.
For formula model objects, it projects EBV and the complete covariance using
all cross-coefficient blocks. Multiple responses produce one result per response.
}
\details{
With P = I - 11'/n, the numerator is EBV'P EBV + trace(P Sigma).
The trace is evaluated as sum(diag(Sigma)) - sum(Sigma)/n; P is not constructed.
The result estimates focal genetic variance conditional on the supplied model
and its components, rather than estimating a population variance component.

EBV and Sigma must refer to the same model, scale and focal quantities.
Identifiers establish alignment, not provenance. Symmetry or positive
semidefiniteness cannot prove that a caller has not filled absent covariances
with zeros. Missing entries must cause an error, not be substituted by zero.

Sigma is tested using Cholesky, with an eigenvalue fallback for semidefinite
matrices. These checks have cubic time complexity. The formula itself has
quadratic time complexity in the number of EBV. Tiny accepted negative
centered traces are clamped to zero and explicitly recorded; substantially
negative values are rejected. This function does not provide a standard error
or confidence interval for VarEst-CF.
}
\value{
A \code{varestcf_result} list. Formula model input returns a \code{table} with
one row per target, named estimates and projection metadata. Vector input returns
\code{estimate}, \code{ebv_term}, \code{error_term}, \code{n},
\code{normalization}, \code{denominator} and \code{diagnostics}. Diagnostics
include alignment, the raw centered trace, roundoff handling, covariance
validation method and elapsed seconds.
}
\seealso{\code{\link{reconstruct_focal_pev}}}
\examples{
ebv <- c(a = 1, b = 2, c = 4)
Sigma <- diag(c(.2, .3, .4))
dimnames(Sigma) <- list(names(ebv), names(ebv))
varest_cf(ebv, Sigma)$estimate
}
