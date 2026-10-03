# Procedencia y referencia matemática / Provenance and mathematical reference

<p align="center">
  <a href="#espanol">Español</a> · <a href="#english">English</a>
</p>

<a id="espanol"></a>

## Español

La implementación pública se organiza alrededor de dos operaciones generales: extraer el bloque focal de la inversa de las ecuaciones completas de un modelo mixto lineal y combinar valores genéticos compatibles con el trazo centrado de su covarianza de error.

El ensamblado y las validaciones están escritos en R usando Matrix. La selección de columnas conserva la incertidumbre de los efectos fijos y de los restantes efectos aleatorios. Los identificadores se utilizan para alinear entradas, sin asumir formatos propios de un software o conjunto de datos.

Los ejemplos y pruebas son sintéticos y autocontenidos. Incluyen contrastes independientes mediante inversión densa y complemento de Schur. La evidencia histórica y los benchmarks con datos privados se mantienen fuera de los archivos públicos.

Las relaciones y los componentes de varianza son entradas explícitas. Los resultados son condicionales a esos componentes; no incluyen su incertidumbre de estimación.

El autor y mantenedor del paquete es Alan Pardo (`pardo.alan@inta.gob.ar`).
El software se distribuye con licencia MIT, de acuerdo con la elección del autor.
La autoría del software y las referencias metodológicas se registran por separado.

### Referencias metodológicas

Sorensen, D., Fernando, R., y Gianola, D. (2001). Inferring the trajectory of genetic
variance in the course of artificial selection. *Genetical Research*, 77(1), 83–94.
[DOI: 10.1017/S0016672300004845](https://doi.org/10.1017/S0016672300004845).

Este trabajo es la referencia metodológica para inferir la varianza genética
de un grupo definido de individuos. El paquete implementa el cálculo con
predicciones y covarianzas completas de error bajo el modelo lineal y los
componentes de varianza suministrados.

Pardo, A. M., Maizón, D. O., Munilla, S., y Legarra, A. (2026). Impact of estimating
genetic variance in the target group on reliability metrics of the Linear
Regression validation method under selection. *Journal of Animal Breeding and
Genetics*, publicación anticipada en línea, 1–12.
[DOI: 10.1111/jbg.70059](https://doi.org/10.1111/jbg.70059).

Este artículo presenta la implementación y evaluación empírica de VarEst-CF
y examina cómo la estimación de la varianza genética focal influye en las
métricas de confiabilidad del método LR. Se publicó en línea el 19 de junio
de 2026. La reconstrucción de PEV/PEC y el estimador del paquete pueden usarse
en modelos mixtos lineales gaussianos sin realizar una validación LR.

`citation("varestcf")` devuelve la cita del software y estas dos referencias
metodológicas. La licencia MIT corresponde al software; los artículos
conservan sus propias condiciones de acceso y reutilización.

---

<a id="english"></a>

## English

The public implementation is organized around two general operations: extracting
the focal block from the inverse of the complete linear mixed-model equations,
and combining compatible breeding values with the centered trace of their
prediction error covariance.

Assembly and validation are written in R using Matrix. Selected columns retain
uncertainty from fixed effects and all remaining random effects. Identifiers
align inputs without assuming formats specific to a software package or dataset.

Examples and tests are synthetic and self-contained. They include independent
comparisons using dense inversion and the Schur complement. Historical evidence
and benchmarks involving private data remain outside the public files.

Relationships and variance components are explicit inputs. Results are
conditional on those components; they do not include estimation uncertainty in
the components.

The package author and maintainer is Alan Pardo (`pardo.alan@inta.gob.ar`).
The software is distributed under the MIT license selected by the author.
Software authorship and methodological references are recorded separately.

### Methodological references

Sorensen, D., Fernando, R., and Gianola, D. (2001). Inferring the trajectory of
genetic variance in the course of artificial selection. *Genetical Research*,
77(1), 83–94.
[DOI: 10.1017/S0016672300004845](https://doi.org/10.1017/S0016672300004845).

This work is the methodological reference for inferring genetic variance in a
defined group of individuals. The package implements the calculation using
predictions and complete prediction error covariances under the supplied linear
model and variance components.

Pardo, A. M., Maizón, D. O., Munilla, S., and Legarra, A. (2026). Impact of
estimating genetic variance in the target group on reliability metrics of the
Linear Regression validation method under selection. *Journal of Animal
Breeding and Genetics*, early online publication, 1–12.
[DOI: 10.1111/jbg.70059](https://doi.org/10.1111/jbg.70059).

This article presents the implementation and empirical evaluation of VarEst-CF
and examines how estimating focal genetic variance affects reliability metrics
of the LR method. It was published online on June 19, 2026. The package's PEV/PEC
reconstruction and estimator can be used with linear Gaussian mixed models
without conducting LR validation.

`citation("varestcf")` returns the software citation and these two methodological
references. The MIT license applies to the software; the articles retain their
own access and reuse conditions.
