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
\emph{English reference below; complete Spanish documentation in the
Español section. Referencia en inglés a continuación; documentación
completa en español en la sección Español.}

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

\section{Español}{
\subsection{Descripción}{
Reconstruye el bloque completo de covarianzas de error de predicción de los
coeficientes aleatorios focales. Forma las ecuaciones completas e identificadas
del modelo mixto, las factoriza una vez y resuelve sólo las columnas necesarias.
Conserva la incertidumbre de efectos fijos y las covarianzas entre coeficientes.
La interfaz de datos y fórmulas prepara diseños, precisiones y EBV compatibles
desde las relaciones y componentes de varianza suministrados.
}
\subsection{Argumentos}{
\describe{
\item{X}{Matriz numérica u objeto Matrix de tipo double con covariables de efectos
fijos identificados. Sus columnas deben ser linealmente independientes. Se
permite un diseño sin columnas.}
\item{Z}{Matriz numérica u objeto Matrix de tipo double con covariables de efectos
aleatorios. Sus nombres de columna identifican coeficientes de forma única. Las
filas deben seguir el mismo orden que X.}
\item{random_precision}{Matriz Q nombrada, simétrica y definida positiva: inversa
de la covarianza conjunta de todos los coeficientes de Z, incluida la escala de
varianza. Filas y columnas se alinean con las columnas de Z. No equivale a una
inversa de relaciones salvo que la escala de varianza sea uno.}
\item{focal_ids}{IDs únicos de animales en la interfaz de fórmulas, o de
coeficientes aleatorios en la interfaz matricial, en el orden de salida deseado.
Los coeficientes pueden representar varios caracteres o regresiones aleatorias.}
\item{residual_precision}{Inversa de la covarianza residual: escalar positivo,
vector positivo por registro o matriz nombrada, simétrica y definida positiva.
Los vectores/matrices nombrados requieren nombres de fila correspondientes en
X/Z y se alinean por ID de registro. Un vector sin nombres se interpreta por
posición.}
\item{block_size}{Entero positivo: número de columnas de la inversa que se
resuelven juntas.}
\item{memory_budget_gib}{Presupuesto orientativo positivo de almacenamiento de
trabajo en GiB. No limita la memoria del sistema operativo ni garantiza RAM
suficiente.}
\item{tolerance}{Tolerancia relativa positiva para simetría y errores hacia atrás
normalizados de las resoluciones por columnas seleccionadas.}
\item{data}{Data frame para la interfaz automática de fórmulas.}
\item{formula}{Respuesta y efectos fijos, por ejemplo \code{weight ~ sex + group}.
Use \code{cbind(y1,y2)} para varias respuestas numéricas. Se admiten respuestas
ausentes.}
\item{animal}{Nombre de la columna de IDs de animales en data.}
\item{relationship}{Matriz de relaciones nombrada y definida positiva: A, H u
otra relación válida para todos los animales del modelo, sin escala de varianza
genética. Aporte exactamente una de relationship o relationship_inverse en la
interfaz de fórmulas con un solo término genético.}
\item{relationship_inverse}{Inversa de relaciones calculada previamente. Se
recomienda para evitar una inversión adicional potencialmente densa. Todos los
IDs de animales deben figurar en filas y columnas. No se procesan pedigrees ni
genotipos.}
\item{random}{Fórmula unilateral de la base de coeficientes aleatorios del
animal; valor predeterminado \code{~1}. La regresión aleatoria, por ejemplo
\code{~1 + time}, admite una sola respuesta.}
\item{genetic_variance}{Varianza genética escalar conocida o matriz de
covarianzas entre coeficientes aleatorios/respuestas. Las matrices nombradas se
alinean por componente.}
\item{residual_variance}{Para una respuesta: escalar positivo conocido, vector
por fila original de data o expresión unilateral evaluada en data. Para varias
respuestas: escalar, vector por respuesta o matriz de covarianzas residuales.}
\item{y}{Vector opcional de respuestas finitas en la interfaz matricial explícita,
en el orden de filas del diseño. Permite calcular EBV focales con las mismas
ecuaciones.}
\item{random_effects}{Lista nombrada alternativa de especificaciones de términos
aleatorios. Cada término requiere \code{id}, una columna de data. Los campos
opcionales son \code{random} (base unilateral, predeterminada \code{~1}),
\code{variance} (escalar o covarianza entre componentes), \code{relationship} o
\code{relationship_inverse}, \code{levels} y \code{allow_missing}. Si no se
aporta una relación, los niveles son independientes y se construye la identidad.
\code{levels} define su orden completo; si se omite, se infieren los niveles
observados. El término genético focal requiere una relación explícitamente
suministrada. No mezcle esta lista con animal, random, genetic_variance ni las
relaciones del nivel superior.}
\item{focal_effect}{Nombre obligatorio del término genético que se resume cuando
se usa random_effects. Todos los demás términos permanecen en las ecuaciones
completas.}
\item{correlated_effects}{Lista opcional de grupos disjuntos, cada uno con
\code{effects} (al menos dos nombres de términos) y \code{covariance} (covarianza
conjunta entre componentes). Los términos del grupo deben compartir la misma
inversa de relaciones y conjunto de IDs. La covarianza sigue el orden de efectos
y luego el de componentes de cada término. Un término escalar usa su nombre;
uno con varios componentes usa término:componente. Las covarianzas nombradas se
alinean por esas etiquetas. Dentro de un grupo puede omitirse la varianza
individual; si se aporta, debe coincidir con su bloque diagonal conjunto. Fuera
de los grupos todas las varianzas son obligatorias.}
}
}
\subsection{Ecuaciones y escala}{
Para diseños X y Z, precisión residual W y precisión aleatoria Q, la matriz C
tiene bloques X'WX, X'WZ, Z'WX y Z'WZ + Q. Si E selecciona los coeficientes
focales dentro del sistema completo, la covarianza es E'C-inversa E. No se
invierte la submatriz principal focal de C ni se consideran conocidos los
efectos fijos.

La covarianza es condicional al modelo y a los componentes suministrados. No
incluye incertidumbre de estimación de componentes. En la interfaz matricial
debe suministrarse una parametrización identificada equivalente de X; las
columnas redundantes se rechazan. En la interfaz de fórmulas se eliminan y
registran columnas fijas redundantes mediante QR dispersa escalada. No se
estiman componentes de varianza.

Con coeficientes agrupados por componente, Q puede construirse mediante
\code{kronecker(solve(G0), Kinv)}. K puede ser A, H u otra covarianza válida.
La interfaz de fórmulas construye Q con la inversa de relaciones y covarianza
genética suministradas; la interfaz matricial recibe Q explícita. Nunca se
sustituye H por A. Las relaciones deben incluir todos los animales del modelo,
no sólo los focales.

Con varias respuestas hay coeficientes fijos separados por carácter y el
apilamiento se realiza por carácter. Con respuestas ausentes y residuos
correlacionados se invierte por separado cada submatriz residual observada.
Los términos aleatorios independientes forman bloques de precisión separados;
para un grupo correlacionado se invierte primero su covarianza conjunta entre
componentes y luego se calcula el producto de Kronecker con la inversa de
relaciones común. Los términos adicionales pueden tener sus propias bases de
regresión aleatoria para una respuesta.

Por defecto se rechazan IDs vacíos o ausentes en filas observadas.
\code{allow_missing = TRUE} asigna explícitamente incidencia cero a ese término;
no infiere ni modela el efecto de un progenitor desconocido. Represente el nivel
desconocido explícitamente cuando su efecto forme parte del modelo previsto.
Se conservan niveles no observados presentes en relaciones o niveles explícitos.

No se automatizan offsets de fórmula, regresión aleatoria multicaracter,
covarianzas generales entre kernels ni bloques residuales multivariados
heterogéneos. El alcance son modelos mixtos lineales gaussianos condicionales a
covarianzas suministradas.
}
\subsection{Recursos y exactitud numérica}{
La referencia valida que Q y una precisión residual general sean definidas
positivas mediante factorizaciones adicionales. Las cotas conservadoras de
ensamblado y los índices CSC de 32 bits pueden rechazar un modelo que una
representación más optimizada admitiría. El llenado del factor y las
asignaciones temporales siguen siendo inciertos; una comprobación posterior a
la factorización no impide todos los casos de falta de memoria. Sigma es densa
y su validación final con Cholesky cuesta tiempo cúbico en el número de
coeficientes focales. No se usan caches ni aproximaciones implícitas.
}
\subsection{Resultado}{
Lista de clase \code{varestcf_fit} con \code{Sigma}, covarianza focal densa y
nombrada; \code{ebv}, valores genéticos compatibles si se aportaron respuestas;
metadatos del modelo de fórmulas cuando corresponde, incluidos el término focal,
los términos del modelo y conteos de incidencias explícitamente nulas;
\code{focal_ids}, su orden; y \code{diagnostics}.

Los diagnósticos incluyen dimensiones, entradas almacenadas, tamaños de objetos,
cota de ensamblado, residuos absolutos y errores hacia atrás por bloque, y
segundos por etapa. Los tamaños y estimaciones no miden el pico de memoria
residente. El error hacia atrás es el máximo residuo absoluto dividido por
la norma infinito de C multiplicada por el máximo valor absoluto de la solución,
más uno. No limita el error hacia adelante en un modelo mal condicionado.
}
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
