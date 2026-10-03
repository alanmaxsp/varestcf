# Capacidad / Computational capacity (reference 0.3)

<p align="center">
  <a href="#espanol">Español</a> · <a href="#english">English</a>
</p>

<a id="espanol"></a>

## Español

### Alcance disponible

- Interfaz automática de datos y fórmulas con K o sus inversas suministradas por el usuario. Prepara diseños, precisión conjunta y EBV; no procesa pedigrees ni genotipos.
- Un carácter con efectos animales o regresiones aleatorias, o varios caracteres con efectos animales correlacionados, incluyendo términos adicionales independientes o correlacionados. Los componentes de varianza deben proporcionarse. Consulte la tabla de cobertura para los límites de preparación automática.
- Si se aporta K sin invertir, la inversión puede ser densa y requiere almacenamiento cuadrático en el número total de animales. Se recomienda suministrar K⁻¹, preferentemente dispersa cuando su estructura lo permita. Una inversa densa sigue siendo costosa.
- Preparar Q = G0⁻¹ ⊗ K⁻¹ puede multiplicar hasta por k² las posiciones de K⁻¹, siendo k el número de coeficientes por animal. Las fórmulas simplifican el uso, pero no eliminan este costo.
- Con varios términos, cada diseño Z ocupa almacenamiento y cada grupo correlacionado tiene su propia precisión conjunta. La cota de preparación suma los diseños, las relaciones y hasta k² veces las posiciones de cada kernel por grupo. Los efectos auxiliares aumentan el número total de ecuaciones aunque no aumenten el bloque focal solicitado.
- Los diseños fijos y la precisión residual se preparan una vez y se comparten entre términos. Aun así, combinar diseños, ensamblar precisiones y permutar bloques puede crear copias temporales.

- Modelos mixtos lineales con diseños explícitos X/Z, precisión residual positiva y precisión conjunta Q definida positiva.
- Precisión conjunta Q suministrada para todos los efectos aleatorios del modelo, siempre con la escala y el orden correctos. El paquete no construye K desde archivos ni procesa genotipos en esta versión.
- Bloque focal completo y cálculo de VarEst-CF. Sin estimación de componentes, caches, GPU, paralelismo distribuido o almacenamiento fuera de RAM.
- Álgebra dispersa proporcionada por Matrix; la referencia es código R transparente sobre ese backend compilado. No incluye optimizaciones C++ propias.

### Costos que deben conocerse

Sean M las ecuaciones completas, r los coeficientes aleatorios, q los coeficientes focales y b el tamaño del bloque de resolución.

| Objeto/operación | Costo relevante |
|---|---|
| Q explícita | Depende de sus posiciones almacenadas; puede ser casi densa |
| C y su factor | Dependen de estructura, ordenación y llenado; no sólo de M |
| B y V | Dos matrices M×b double: aproximadamente 16Mb bytes, más temporales |
| Sigma | q×q double: 8q² bytes, más copias/validaciones |
| Validación de Q y W matricial | Factorizaciones adicionales en esta referencia estricta |
| Validación final de Sigma | Cholesky denso, O(q³) |
| Fórmula VarEst-CF | O(q²) para q EBV; la validación de la matriz añade O(q³) |

Una Sigma de 10.000 coeficientes requiere 0,745 GiB sólo en valores; de 20.000, 2,980 GiB. Estos valores no son presupuestos recomendados. La función estadística recibe una matriz ya construida y no impone un presupuesto propio; el usuario debe poder almacenarla y validarla.

### Guardas implementadas

`memory_budget_gib` comprueba estimaciones de salida/buffers, inputs, validaciones, ensamblado y almacenamiento tras factorizar C. Las constantes de trabajo incluyen varias copias previstas, pero no garantizan capturar todos los temporales de R/CHOLMOD. Los objetos de entrada ya existen antes de entrar en la función; el presupuesto no controla su creación ni la RAM ocupada por otras aplicaciones.

El ensamblado R puede usar matrices CSC generales, con hasta dos triángulos almacenados temporalmente. Para precisión residual diagonal, se usa una cota estructural a partir de las columnas activas por fila del diseño. Para W general se adopta la cota conservadora M². Una cota excesiva puede rechazar un problema que otro backend admitiría. No debe aumentarse el presupuesto sin evaluar recursos.

La representación actual usa índices CSC de 32 bits: se rechazan cotas de almacenamiento de 2.147.483.647 posiciones o más, y dimensiones fuera del rango admitido. Este límite depende de la representación, no de la RAM disponible. No equivale a un máximo fijo de animales.

El llenado del factor no se conoce con certeza antes de calcularlo. Las guardas son orientativas y portables porque usan dimensiones/objetos, **no una detección garantizada de RAM disponible**. No se promete impedir todas las terminaciones del sistema por falta de memoria. Los residuos normalizados son errores hacia atrás; no garantizan error hacia adelante pequeño si el sistema está mal condicionado.

### Capacidad demostrada y benchmarks pendientes

La suite de aceptación utiliza ocho animales sintéticos con relaciones genealógicas, una corrección genómica inventada, registros repetidos y hasta cuatro coeficientes por animal. Sirve para comprobar exactitud y contratos, no para certificar escalabilidad.

Antes de anunciar capacidades públicas se requiere un benchmark sintético escalonado que registre por separado construcción de Q, validaciones, ensamblado, factorización, extracción y estadístico. Debe variar el número de animales, la densidad de la precisión, los componentes y términos aleatorios, los efectos fijos y los focales, y registrar hardware, hilos, versiones y pico RSS medido externamente. Una única matriz grande no establece una regla universal.

La ruta optimizada futura deberá reproducir esta referencia con tolerancias justificadas. Los límites de recursos no autorizan a modificar las relaciones, omitir efectos presentes en el modelo ni reemplazar covarianzas ausentes por cero.

---

<a id="english"></a>

## English

### Available scope

- Automatic data/formula interface with user-supplied K or its inverse. It
  prepares designs, joint precision and EBV; it does not process pedigrees or
  genotypes.
- One trait with animal effects or random regressions, or multiple traits with
  correlated animal effects, including additional independent or correlated
  terms. Variance components must be supplied. See the
  [coverage table](MODELOS.md#english) for limits of automatic preparation.
- If K is supplied without inversion, its inverse may be dense and require
  storage quadratic in the total number of animals. Supplying K⁻¹ is recommended,
  preferably sparse when its structure permits. A dense inverse remains costly.
- Preparing Q = G0⁻¹ ⊗ K⁻¹ may multiply the stored entries of K⁻¹ by up to k²,
  where k is the number of coefficients per animal. Formulas simplify usage but
  do not remove this cost.
- With multiple terms, each Z design uses storage and each correlated group has
  its own joint precision. The preparation bound sums designs, relationships and
  up to k² times the entries of each kernel per group. Auxiliary effects increase
  the total number of equations even when they do not increase the requested
  focal block.
- Fixed designs and residual precision are prepared once and shared across
  terms. Combining designs, assembling precision matrices and permuting blocks
  may still create temporary copies.
- Linear mixed models with explicit X/Z designs, positive residual precision and
  positive definite joint random precision Q.
- Joint Q supplied for all random effects in the model, with the correct scale
  and order. This version does not build K from files or process genotypes.
- Complete focal block and VarEst-CF calculation. No component estimation,
  caches, GPU, distributed parallelism or storage outside RAM.
- Sparse algebra is provided by Matrix; the reference is transparent R code
  using that compiled backend. It includes no custom C++ optimizations.

### Computational costs

Let M be the number of complete equations, r the number of random coefficients,
q the number of focal coefficients, and b the solve block size.

| Object/operation | Relevant cost |
|---|---|
| Explicit Q | Depends on its stored entries; it may be nearly dense |
| C and its factor | Depend on structure, ordering and fill, not only M |
| B and V | Two M×b double matrices: approximately 16Mb bytes, plus temporaries |
| Sigma | q×q double matrix: 8q² bytes, plus copies/validation |
| Validation of Q and matrix W | Additional factorizations in this strict reference |
| Final validation of Sigma | Dense Cholesky, O(q³) |
| VarEst-CF formula | O(q²) for q EBV; matrix validation adds O(q³) |

A Sigma with 10,000 coefficients needs 0.745 GiB for its values alone; with
20,000, 2.980 GiB. These figures are not recommended budgets. The statistical
function receives an existing matrix and does not impose its own budget; users
must be able to store and validate it.

### Implemented guards

`memory_budget_gib` checks estimates for output/buffers, inputs, validation,
assembly and storage after factoring C. Working constants allow for several
expected copies, but cannot guarantee coverage of every R/CHOLMOD temporary.
Input objects already exist before the function starts; the budget does not
control their creation or the RAM used by other applications.

R assembly may use general CSC matrices, with both triangles temporarily stored.
For diagonal residual precision, a structural bound uses the active design
columns per row. For general W, the conservative M² bound is used. An excessive
bound may reject a problem that another backend could handle. Do not increase
the budget without assessing resources.

The current representation uses 32-bit CSC indices: storage bounds of
2,147,483,647 entries or more, and dimensions outside the supported range, are
rejected. This limit depends on representation rather than available RAM. It
does not imply a fixed maximum number of animals.

Factor fill is not known with certainty before factorization. Guards are
advisory and portable because they use dimensions/object sizes,
**not guaranteed detection of available RAM**. They do not promise to prevent
all out-of-memory termination. Normalized residuals are backward errors; they
do not guarantee small forward error for an ill-conditioned system.

### Demonstrated capacity and pending benchmarks

The acceptance suite uses eight synthetic animals with pedigree relationships,
an invented genomic correction, repeated records and up to four coefficients
per animal. It checks numerical accuracy and input contracts, not scalability.

Before announcing public capacity claims, a staged synthetic benchmark must
record Q construction, validation, assembly, factorization, extraction and the
statistic separately. It must vary animal numbers, precision density, random
components and terms, fixed effects and focal animals, and record hardware,
threads, versions and externally measured peak RSS. A single large matrix does
not establish a universal rule.

A future optimized path must reproduce this reference with justified tolerances.
Resource limits do not justify changing relationships, omitting model effects
or replacing absent covariances with zeros.
