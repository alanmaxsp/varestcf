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
\emph{English reference below; complete Spanish documentation in the
Español section. Referencia en inglés a continuación; documentación
completa en español en la sección Español.}

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

\section{Español}{
\subsection{Descripción}{
Calcula la suma de cuadrados de EBV centrados y el trazo centrado de la covarianza
de error de predicción, dividida por el denominador elegido. Esta función no
reconstruye ecuaciones de modelos mixtos ni pedigrees. Para objetos de modelos
de fórmulas proyecta EBV y la covarianza completa utilizando todos los bloques
cruzados entre coeficientes. Con varias respuestas devuelve un resultado por
respuesta.
}
\subsection{Argumentos}{
\describe{
\item{ebv}{Objeto de modelo reconstruido con EBV, o vector numérico finito y
nombrado de al menos dos valores genéticos focales.}
\item{Sigma}{Matriz densa explícitamente completa de covarianzas de error de
predicción, con IDs únicos en filas y columnas correspondientes a los nombres
de EBV. Las entradas se alinean por ID. Se rechazan entradas Matrix dispersas
para evitar interpretar omisiones estructurales como covarianzas cero
suministradas. Omita Sigma cuando aporta un objeto de modelo.}
\item{normalization}{Denominador n (predeterminado) o n menos uno.}
\item{tolerance}{Tolerancia relativa positiva para validar la covarianza y el
redondeo numérico del trazo centrado.}
\item{at}{Puntos de evaluación de regresión aleatoria: vector numérico para una
variable o data frame cuyas columnas correspondan a la fórmula aleatoria.}
\item{component}{Nombres opcionales de coeficientes individuales/respuestas a
evaluar. Se excluye mutuamente con at. La regresión aleatoria requiere at o
component.}
}
}
\subsection{Ecuación y condiciones}{
Con P = I - 11'/n, el numerador es EBV'P EBV + trace(P Sigma). El trazo se evalúa
como sum(diag(Sigma)) - sum(Sigma)/n; P no se construye. Se divide por n o n-1.
El resultado estima la varianza genética focal condicional al modelo y sus
componentes suministrados; no estima un componente de varianza poblacional.

EBV y Sigma deben corresponder al mismo modelo, escala y cantidades focales.
Los IDs establecen alineación, no procedencia. Simetría o semidefinición positiva
no demuestran que el usuario no haya completado con ceros covarianzas ausentes.
Las entradas faltantes deben causar un error, no sustituirse por cero.

Sigma se valida con Cholesky y, para matrices semidefinidas, con autovalores
como alternativa. Estas comprobaciones cuestan tiempo cúbico. La fórmula cuesta
tiempo cuadrático en el número de EBV. Trazos centrados negativos pequeños,
aceptados por tolerancia, se ajustan a cero y se registran explícitamente;
valores sustancialmente negativos se rechazan. Esta función no proporciona
error estándar ni intervalo de confianza de VarEst-CF.
}
\subsection{Resultado}{
Lista de clase \code{varestcf_result}. Para un modelo de fórmulas, devuelve
\code{table} con una fila por objetivo, estimaciones nombradas y metadatos de
proyección. Para un vector, devuelve \code{estimate}, \code{ebv_term},
\code{error_term}, \code{n}, \code{normalization}, \code{denominator} y
\code{diagnostics}. Los diagnósticos incluyen alineación, trazo centrado sin
ajuste, tratamiento del redondeo, método de validación de covarianza y segundos
transcurridos.
}
}

\seealso{\code{\link{reconstruct_focal_pev}}}
\examples{
ebv <- c(a = 1, b = 2, c = 4)
Sigma <- diag(c(.2, .3, .4))
dimnames(Sigma) <- list(names(ebv), names(ebv))
varest_cf(ebv, Sigma)$estimate
}
