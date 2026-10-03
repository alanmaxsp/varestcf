# VarEst-CF

<p align="center">
  <a href="#espanol">Español</a> · <a href="#english">English</a>
</p>

<a id="espanol"></a>

## Español

Reconstrucción del bloque completo de error de predicción de valores genéticos focales y cálculo de su varianza mediante VarEst-CF, para modelos mixtos lineales.

Dos funciones públicas:

- `reconstruct_focal_pev()` prepara las ecuaciones, calcula predicciones compatibles y obtiene todas las PEV/PEC focales.
- `varest_cf()` calcula el estadístico a partir de esas predicciones y su matriz completa de error.

### Instalación

Recomendamos usar una versión actual de R. El paquete declara R >= 4.1.0 y
utiliza Matrix; `remotes` se usa sólo para instalar desde GitHub.

Una vez publicado el repositorio, ejecute en la consola de R:

```r
install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github(
  "alanmaxsp/varestcf",
  upgrade = "never",
  build = FALSE,
  repos = "https://cloud.r-project.org"
)
library(varestcf)
```

Las dependencias necesarias se instalan si faltan. Este comando instala la rama
principal; `upgrade = "never"` evita actualizar otras dependencias por iniciativa
del instalador. El argumento `build = FALSE` permite instalar
directamente la fuente de este paquete R sin construir documentación adicional.

### Versión estable, desarrollo y versiones numeradas

La instalación anterior usa `main`, la rama reservada para cambios revisados.
`develop` es la rama para cambios y pruebas de la próxima versión. Sus cambios
no modifican `main` hasta que se revisan y se incorporan explícitamente.

Para probar la versión de desarrollo, elíjala expresamente:

```r
remotes::install_github("alanmaxsp/varestcf@develop", build = FALSE)
```

Para reproducir un análisis con una versión numerada, use su etiqueta. Este
ejemplo será válido **cuando publiquemos la etiqueta `v0.3.0`**:

```r
remotes::install_github("alanmaxsp/varestcf@v0.3.0", build = FALSE)
```

Todavía no hay versiones numeradas publicadas. Una rama puede avanzar; una
etiqueta publicada se conserva sin cambios y las correcciones reciben una nueva
versión. Registre `packageVersion("varestcf")` y `sessionInfo()`; si usa una rama,
registre además el commit instalado (`packageDescription("varestcf")$RemoteSha`).

En un repositorio público, `develop` también es visible y descargable, pero sólo
se instala si el usuario la selecciona. Instalar estable y desarrollo en la misma
biblioteca de R reemplaza la instalación anterior; use bibliotecas separadas si
necesita conservar ambas simultáneamente.

La [guía de trabajo](docs/DESARROLLO.md#espanol) explica cómo cambiar de rama en GitHub Desktop
y revisar cambios antes de incorporarlos a `main`.

### Primer análisis completo

Copie este ejemplo después de instalar el paquete. Todos los datos son
inventados. K es una matriz de relaciones suministrada; el paquete no necesita
construirla desde un pedigree.

```r
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
```

El resultado esperado, salvo redondeo numérico, es:

```text
VarEst-CF
      target n  estimate  ebv_term error_term normalization
 (Intercept) 3 0.9779592 0.6160544  0.3619048             n
```

`estimate` es la varianza genética del grupo focal estimada mediante VarEst-CF.
`ebv_term` es la contribución de sus EBV centrados y `error_term` es la corrección
obtenida de la matriz completa de PEV/PEC. La normalización predeterminada
divide por el número de animales focales.

El grupo focal incluye `animal4`, que no tiene observaciones. Se conserva porque
está incluido en K y en las ecuaciones. `fit$Sigma` es una matriz completa 3×3:
su diagonal contiene las PEV y las entradas fuera de la diagonal contienen las
PEC entre los focales. `fit$ebv` sigue el mismo orden.

El mismo ejemplo está instalado con el paquete y se puede ejecutar mediante:

```r
source(system.file("examples", "getting_started.R", package = "varestcf"))
```

### Con tus datos

El usuario aporta datos, fórmulas, componentes de varianza y una matriz de relaciones K o su inversa. K puede ser genealógica (A), genómica (G), combinada (H) u otra matriz válida para el efecto genético especificado.

En las plantillas siguientes, sustituya `records`, `Kinv`, `focal_ids` y los
componentes de varianza por las entradas de su modelo. Las columnas nombradas
en las fórmulas deben existir en `records`.

```r
library(varestcf)

fit <- reconstruct_focal_pev(
  data = records,
  formula = trait ~ sex + age + group,
  animal = "animal",
  relationship_inverse = Kinv,
  focal_ids = focal_ids,
  genetic_variance = 4,
  residual_variance = 5
)
varest_cf(fit)
```

Los valores 4 y 5 son ilustrativos. Deben corresponder al modelo cuyos valores genéticos se quieren analizar.

La matriz de relaciones debe tener IDs únicos en filas y columnas e incluir todos los animales del modelo, también los que no tienen observaciones. La función alinea los identificadores. Aporte exactamente una entrada: `relationship_inverse = Kinv` o `relationship = K`. Se recomienda la inversa para evitar una inversión adicional potencialmente densa.

La escala genética se aporta por separado: K es la matriz de relaciones; `genetic_variance` es la varianza o covarianza entre componentes genéticos. La función identifica y registra columnas fijas redundantes, conserva la incertidumbre de los efectos fijos y obtiene predicciones del mismo modelo que la matriz focal.

Consulte `fit$Sigma`, `fit$ebv` y `fit$diagnostics` para acceder a la matriz completa, las predicciones y las comprobaciones.

### Regresión aleatoria

```r
fit_rr <- reconstruct_focal_pev(
  data = records, formula = trait ~ sex + group,
  animal = "animal", relationship_inverse = Kinv,
  focal_ids = focal_ids,
  random = ~ 1 + time,
  genetic_variance = G0,
  residual_variance = records$residual_var
)
varest_cf(fit_rr, at = c(10, 20, 30))
```

G0 es la covarianza entre coeficientes de la base aleatoria. Los puntos se expresan en la misma escala usada en la fórmula. Las proyecciones incluyen todos los cruces entre coeficientes y animales.

La varianza residual puede ser constante, un vector por fila original o una expresión conocida, por ejemplo `residual_variance = ~ exp(b0 + b1 * time)`. Para clases, prepare el vector mediante la clase de cada observación. La función recibe las varianzas finales, no estima sus parámetros. La heterogeneidad modifica las ecuaciones y debe reproducir la estructura del modelo completo.

### Múltiples caracteres y otros modelos

Para varios caracteres use `formula = cbind(trait1, trait2) ~ sex + group`, una matriz genética entre caracteres y un vector de varianzas residuales o una matriz residual. Las respuestas ausentes pueden ser `NA`; los predictores utilizados deben estar completos. `varest_cf(fit)` devuelve un resultado por carácter.

La interfaz admite varios términos aleatorios. Todos permanecen en las ecuaciones; se selecciona explícitamente el efecto genético focal. Omitir otros términos presentes en el modelo cambia sus PEV/PEC.

### Varios términos aleatorios

Para un modelo de repetibilidad:

```r
fit <- reconstruct_focal_pev(
  data = records, formula = trait ~ sex + group,
  random_effects = list(
    additive = list(id = "animal", relationship_inverse = Kinv, variance = 4),
    permanent = list(id = "animal", variance = 2)
  ),
  focal_effect = "additive",
  focal_ids = focal_ids,
  residual_variance = 5
)
varest_cf(fit)
```

Sin una matriz de relaciones, los niveles del término se tratan como independientes. Para camada u otro efecto común basta con usar la columna de identificación correspondiente. Cada término también puede tener una base `random = ~ 1 + time` y su matriz de covarianzas para una respuesta.

Para efectos genéticos directos y maternos correlacionados:

```r
fit_dm <- reconstruct_focal_pev(
  data = records, formula = trait ~ sex + group,
  random_effects = list(
    direct = list(id = "animal", relationship_inverse = Kinv),
    maternal = list(id = "dam", relationship_inverse = Kinv)
  ),
  correlated_effects = list(
    list(effects = c("direct", "maternal"), covariance = Gdm)
  ),
  focal_effect = "direct",
  focal_ids = focal_ids,
  residual_variance = 5
)
varest_cf(fit_dm)
```

Gdm es la matriz conjunta 2×2 de varianzas y covarianza directa/materna. Los efectos correlacionados deben compartir la misma matriz de relaciones, con los mismos IDs. La función invierte la covarianza conjunta completa.

Las IDs de los términos deben estar completas en las filas observadas. La opción avanzada `allow_missing = TRUE` significa explícitamente incidencia cero para ese término; no representa automáticamente un efecto parental desconocido.

La [tabla de cobertura](docs/MODELOS.md#espanol) distingue capacidades implementadas, extensiones que requieren matrices explícitas y modelos fuera del alcance actual.

### Otros ejemplos reproducibles

Después de instalar el paquete, puede ejecutar otros ejemplos completos:

```r
source(system.file("examples", "formula_workflow.R", package = "varestcf"))
source(system.file("examples", "multiple_effects.R", package = "varestcf"))
```

Los ejemplos son completamente sintéticos. Dependencias: R y Matrix.

La [guía de desarrollo](docs/DESARROLLO.md#espanol) explica la suite consolidada, los
checks en Windows/Linux y la revisión de archivos antes de publicar.

### Recursos y referencia matemática

La capacidad depende de la densidad de la precisión, el llenado del factor y el número de coeficientes focales. El presupuesto orientativo predeterminado es 1 GiB, configurable con `memory_budget_gib`. Consulte [capacidad y límites](docs/CAPACIDAD.md#espanol) para dimensionar el análisis.

La interfaz avanzada recibe X, Z, la precisión conjunta de todos los efectos aleatorios y la precisión residual. También puede calcularse `varest_cf(ebv, Sigma)` con valores genéticos y una matriz completa compatible. La ayuda de las funciones documenta ecuaciones, escalas y validaciones en ambos idiomas:

```r
help("reconstruct_focal_pev", package = "varestcf")
help("varest_cf", package = "varestcf")
```

El alcance es modelos mixtos lineales gaussianos, con componentes suministrados.

Versión 0.3.0 experimental. Autor y mantenedor: Alan Pardo
(`pardo.alan@inta.gob.ar`). Licencia [MIT](LICENSE.md).

`citation("varestcf")` devuelve la cita del software y las referencias
metodológicas: [Sorensen et al. (2001)](https://doi.org/10.1017/S0016672300004845)
y [Pardo et al. (2026)](https://doi.org/10.1111/jbg.70059).
Su relación con el paquete se explica en [procedencia](docs/PROVENIENCIA.md#espanol).

---

<a id="english"></a>

## English

Reconstruct the complete prediction error covariance block of focal breeding
values and estimate their genetic variance using VarEst-CF for linear mixed models.

Two public functions:

- `reconstruct_focal_pev()` prepares the equations, computes compatible predictions
  and obtains all focal PEV/PEC entries.
- `varest_cf()` computes the statistic from those predictions and their complete
  prediction error covariance matrix.

### Installation

We recommend a current version of R. The package declares R >= 4.1.0 and uses
Matrix; `remotes` is used only to install from GitHub.

Once the repository is public, run the following in the R console:

```r
install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github(
  "alanmaxsp/varestcf",
  upgrade = "never",
  build = FALSE,
  repos = "https://cloud.r-project.org"
)
library(varestcf)
```

Required dependencies are installed if missing. This command installs the main
branch; `upgrade = "never"` prevents the installer from updating other dependencies
on its own initiative. The `build = FALSE` argument installs this R package
directly from its source without building additional documentation.

### Stable, development and numbered versions

The installation above uses `main`, the branch reserved for reviewed changes.
`develop` is the branch for changes and tests for the next version. Its changes
do not modify `main` until they are reviewed and explicitly merged.

To try the development version, select it explicitly:

```r
remotes::install_github("alanmaxsp/varestcf@develop", build = FALSE)
```

To reproduce an analysis with a numbered version, use its tag. This example will
work **once we publish the `v0.3.0` tag**:

```r
remotes::install_github("alanmaxsp/varestcf@v0.3.0", build = FALSE)
```

No numbered versions have been published yet. A branch may advance; a published
tag is retained unchanged and fixes receive a new version. Record
`packageVersion("varestcf")` and `sessionInfo()`; when using a branch, also record the
installed commit (`packageDescription("varestcf")$RemoteSha`).

In a public repository, `develop` is also visible and downloadable, but is only
installed when users select it. Installing stable and development in the same R
library replaces the previous installation; use separate libraries if you need
both simultaneously.

The [working guide](docs/DESARROLLO.md#english) explains switching branches in GitHub Desktop
and reviewing changes before merging them into `main`.

### First complete analysis

Copy this example after installing the package. All data are synthetic. K is a
supplied relationship matrix; the package does not need to construct it from a
pedigree.

```r
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
```

The expected result, apart from numerical rounding, is:

```text
VarEst-CF
      target n  estimate  ebv_term error_term normalization
 (Intercept) 3 0.9779592 0.6160544  0.3619048             n
```

`estimate` is the focal group's genetic variance estimated by VarEst-CF.
`ebv_term` is the contribution of centered EBV, and `error_term` is the correction
obtained from the complete PEV/PEC matrix. The default normalization divides by
the number of focal animals.

The focal group includes `animal4`, which has no observations. It is retained
because it is included in K and in the equations. `fit$Sigma` is a complete 3×3
matrix: its diagonal contains PEV and its off-diagonal entries contain PEC between
focal animals. `fit$ebv` follows the same order.

The same example is installed with the package and can be run with:

```r
source(system.file("examples", "getting_started.R", package = "varestcf"))
```

### With your data

The user supplies data, formulas, variance components and a relationship matrix K
or its inverse. K may be a pedigree (A), genomic (G), combined (H), or another valid
relationship matrix for the specified genetic effect.

In the following templates, replace `records`, `Kinv`, `focal_ids` and the variance
components with the inputs for your model. Columns named in the formulas must
exist in `records`.

```r
library(varestcf)

fit <- reconstruct_focal_pev(
  data = records,
  formula = trait ~ sex + age + group,
  animal = "animal",
  relationship_inverse = Kinv,
  focal_ids = focal_ids,
  genetic_variance = 4,
  residual_variance = 5
)
varest_cf(fit)
```

The values 4 and 5 are illustrative. They must match the model whose breeding
values are being analyzed.

The relationship matrix must have unique row and column IDs and include all
animals in the model, including those without observations. The function aligns
identifiers. Supply exactly one input: `relationship_inverse = Kinv` or
`relationship = K`. The inverse is preferred to avoid an additional, potentially
dense inversion.

Genetic scale is supplied separately: K is the relationship matrix;
`genetic_variance` is the variance or covariance between genetic components.
The function identifies and reports redundant fixed-effect columns, retains
fixed-effect uncertainty, and obtains predictions from the same model as the
focal covariance matrix.

Use `fit$Sigma`, `fit$ebv` and `fit$diagnostics` to access the complete matrix,
predictions and numerical checks.

### Random regression

```r
fit_rr <- reconstruct_focal_pev(
  data = records, formula = trait ~ sex + group,
  animal = "animal", relationship_inverse = Kinv,
  focal_ids = focal_ids,
  random = ~ 1 + time,
  genetic_variance = G0,
  residual_variance = records$residual_var
)
varest_cf(fit_rr, at = c(10, 20, 30))
```

G0 is the covariance between coefficients of the random basis. Evaluation points
use the same scale as the formula. Projections include all cross-covariances
between coefficients and animals.

Residual variance can be constant, a vector per original data row, or a known
expression, for example `residual_variance = ~ exp(b0 + b1 * time)`. For classes,
prepare the vector using each observation's class. The function receives the
final variances; it does not estimate their parameters. Heterogeneity changes
the equations and must reproduce the complete model's structure.

### Multiple traits and other models

For multiple traits, use `formula = cbind(trait1, trait2) ~ sex + group`, a genetic
covariance matrix between traits, and a vector of residual variances or a residual
covariance matrix. Missing responses may be `NA`; predictors used in the model
must be complete. `varest_cf(fit)` returns one result per trait.

The interface supports multiple random terms. All remain in the equations; the
focal genetic effect is selected explicitly. Omitting other terms present in the
model changes its PEV/PEC.

### Multiple random terms

For a repeatability model:

```r
fit <- reconstruct_focal_pev(
  data = records, formula = trait ~ sex + group,
  random_effects = list(
    additive = list(id = "animal", relationship_inverse = Kinv, variance = 4),
    permanent = list(id = "animal", variance = 2)
  ),
  focal_effect = "additive",
  focal_ids = focal_ids,
  residual_variance = 5
)
varest_cf(fit)
```

Without a relationship matrix, a term's levels are treated as independent. For
litter or another common effect, use the corresponding identification column.
Each term may also have a `random = ~ 1 + time` basis and its coefficient
covariance matrix for a single response.

For correlated direct and maternal genetic effects:

```r
fit_dm <- reconstruct_focal_pev(
  data = records, formula = trait ~ sex + group,
  random_effects = list(
    direct = list(id = "animal", relationship_inverse = Kinv),
    maternal = list(id = "dam", relationship_inverse = Kinv)
  ),
  correlated_effects = list(
    list(effects = c("direct", "maternal"), covariance = Gdm)
  ),
  focal_effect = "direct",
  focal_ids = focal_ids,
  residual_variance = 5
)
varest_cf(fit_dm)
```

Gdm is the joint 2×2 matrix of direct and maternal variances and their covariance.
Correlated effects must share the same relationship matrix, with the same IDs.
The function inverts the complete joint covariance.

Term IDs must be complete on observed rows. The advanced option
`allow_missing = TRUE` explicitly means zero incidence for that term; it does
not automatically represent an unknown parental effect.

The [model coverage table](docs/MODELOS.md#english) distinguishes implemented
capabilities, extensions requiring explicit matrices, and models outside the
current scope.

### More reproducible examples

After installing the package, you can run other complete examples:

```r
source(system.file("examples", "formula_workflow.R", package = "varestcf"))
source(system.file("examples", "multiple_effects.R", package = "varestcf"))
```

All examples are synthetic. Dependencies: R and Matrix.

The [development guide](docs/DESARROLLO.md#english) explains the consolidated
suite, Windows/Linux checks and file review before publication.

### Resources and mathematical reference

Capacity depends on precision density, factor fill and the number of focal
coefficients. The default advisory budget is 1 GiB, configurable through
`memory_budget_gib`. See [capacity and limits](docs/CAPACIDAD.md#english) to plan
your analysis.

The advanced interface receives X, Z, the joint precision of all random effects,
and residual precision. You can also use `varest_cf(ebv, Sigma)` with breeding
values and a compatible complete covariance matrix. Function help documents
equations, scales and validations in both languages:

```r
help("reconstruct_focal_pev", package = "varestcf")
help("varest_cf", package = "varestcf")
```

The scope is linear Gaussian mixed models with supplied variance components.

Experimental version 0.3.0. Author and maintainer: Alan Pardo
(`pardo.alan@inta.gob.ar`). License: [MIT](LICENSE.md).

`citation("varestcf")` returns the software citation and methodological references:
[Sorensen et al. (2001)](https://doi.org/10.1017/S0016672300004845) and
[Pardo et al. (2026)](https://doi.org/10.1111/jbg.70059).
Their relationship to the package is explained in
[provenance](docs/PROVENIENCIA.md#english).
