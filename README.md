# VarEst-CF

Reconstrucción del bloque completo de error de predicción de valores genéticos focales y cálculo de su varianza mediante VarEst-CF, para modelos mixtos lineales.

Dos funciones públicas:

- `reconstruct_focal_pev()` prepara las ecuaciones, calcula predicciones compatibles y obtiene todas las PEV/PEC focales.
- `varest_cf()` calcula el estadístico a partir de esas predicciones y su matriz completa de error.

## Uso habitual

El usuario aporta datos, fórmulas, componentes de varianza y una matriz de relaciones K o su inversa. K puede ser genealógica (A), genómica (G), combinada (H) u otra matriz válida para el efecto genético especificado.

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

## Regresión aleatoria

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

## Múltiples caracteres y otros modelos

Para varios caracteres use `formula = cbind(trait1, trait2) ~ sex + group`, una matriz genética entre caracteres y un vector de varianzas residuales o una matriz residual. Las respuestas ausentes pueden ser `NA`; los predictores utilizados deben estar completos. `varest_cf(fit)` devuelve un resultado por carácter.

La interfaz admite varios términos aleatorios. Todos permanecen en las ecuaciones; se selecciona explícitamente el efecto genético focal. Omitir otros términos presentes en el modelo cambia sus PEV/PEC.

## Varios términos aleatorios

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

La [tabla de cobertura](docs/MODELOS.md) distingue capacidades implementadas, extensiones que requieren matrices explícitas y modelos fuera del alcance actual.

## Ejemplo reproducible

Instale el paquete local con `R CMD INSTALL varestcf_0.3.0.tar.gz` y ejecute:

```r
source(system.file("examples", "formula_workflow.R", package = "varestcf"))
source(system.file("examples", "multiple_effects.R", package = "varestcf"))
```

Los ejemplos son completamente sintéticos. Dependencias: R y Matrix.

La [guía de desarrollo](docs/DESARROLLO.md) explica la suite consolidada, los
checks en Windows/Linux y la revisión de archivos antes de publicar.

## Recursos y referencia matemática

La capacidad depende de la densidad de la precisión, el llenado del factor y el número de coeficientes focales. El presupuesto orientativo predeterminado es 1 GiB, configurable con `memory_budget_gib`. Consulte [capacidad y límites](docs/CAPACIDAD.md) para dimensionar el análisis.

La interfaz avanzada recibe X, Z, la precisión conjunta de todos los efectos aleatorios y la precisión residual. También puede calcularse `varest_cf(ebv, Sigma)` con valores genéticos y una matriz completa compatible. La ayuda de las funciones documenta ecuaciones, escalas y validaciones.

El alcance es modelos mixtos lineales gaussianos, con componentes suministrados.

Versión 0.3.0 experimental. Autor y mantenedor: Alan Pardo
(`pardo.alan@inta.gob.ar`). Licencia [MIT](LICENSE.md).

`citation("varestcf")` devuelve la cita del software y las referencias
metodológicas: [Sorensen et al. (2001)](https://doi.org/10.1017/S0016672300004845)
y [Pardo et al. (2026)](https://doi.org/10.1111/jbg.70059).
Su relación con el paquete se explica en [procedencia](docs/PROVENIENCIA.md).
