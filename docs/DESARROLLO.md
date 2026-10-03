# Desarrollo y comprobaciones / Development and checks

<p align="center">
  <a href="#espanol">Español</a> · <a href="#english">English</a>
</p>

<a id="espanol"></a>

## Español

### Contenido del repositorio

El repositorio contiene la referencia R, ayuda de las dos funciones públicas,
documentación, ejemplos sintéticos, una suite consolidada de aceptación,
metadatos de autoría/licencia/citas y configuración de GitHub Actions.

Las relaciones, observaciones e identificadores de los ejemplos se generan en
`inst/examples/synthetic_models.R`. La inversión densa de ese archivo es un
oráculo para problemas pequeños, separado del algoritmo de extracción por bloques
del paquete. No es un método recomendado para modelos grandes.

Los datos reales, los benchmarks privados, los diagnósticos históricos y los
archivos generados se conservan fuera del repositorio. No forman parte de la
instalación ni son necesarios para ejecutar las pruebas públicas.

El [primer análisis del README](../README.md#primer-análisis-completo) coincide
con `inst/examples/getting_started.R`: se puede copiar completo y también ejecutar
desde la instalación del paquete. Su matriz de relaciones y sus observaciones
son explícitas; no depende de archivos auxiliares de los otros ejemplos.

### Chequeo automático

El flujo [R-CMD-check](../.github/workflows/R-CMD-check.yaml) se ejecuta en cada
push, pull request y ejecución manual, con R release en Windows y Ubuntu.

1. Revisa rutas, formatos, tamaño y contenido binario de los archivos registrados
   en Git mediante [el control de archivos](../.github/scripts/check-public-files.R).
2. Instala las dependencias del paquete y `rcmdcheck`.
3. Construye e instala el paquete y ejecuta `R CMD check --no-manual`, incluyendo
   la ayuda, los ejemplos de ayuda y `tests/acceptance.R`.
4. Guarda los resultados y muestra el final de la suite de aceptación.

Los errores y las advertencias hacen fallar el chequeo. Las notas deben revisarse.
La suite numérica se ejecuta una sola vez por sistema, dentro de R CMD check:
no se vuelve a ejecutar como paso adicional.

Las acciones externas están fijadas por commit y el flujo usa permisos de lectura.
Sus parámetros se basan en las instrucciones oficiales de
[r-lib/actions](https://github.com/r-lib/actions/tree/v2/check-r-package) y
[actions/checkout](https://github.com/actions/checkout). Al actualizar una acción,
revise su versión y cambie el commit fijado en el flujo.

Estos trabajos contrastan portabilidad con las versiones que instalan en cada
ejecución. No fijan una única versión de R/Matrix para todos los análisis.
Los informes guardan las versiones efectivamente usadas.
La compatibilidad con todas las versiones de R desde el mínimo declarado
no queda certificada por esta matriz de dos sistemas.

### Revisión de publicación

Desde la raíz del repositorio, con las fuentes registradas o preparadas en Git:

```text
Rscript .github/scripts/check-public-files.R
```

El control admite únicamente los formatos de texto usados por esta versión,
rechaza formatos no revisados y archivos mayores de 1 MiB, y bloquea directorios
de datos, outputs, benchmarks y caches. La lista es deliberadamente pequeña:
una incorporación que requiera otro formato debe revisarse antes de ampliarla.

El control no puede determinar si una matriz escrita dentro de un script procede
de datos reales ni si existen derechos de redistribución. La revisión de
procedencia sigue siendo necesaria. `.gitignore` evita incorporaciones
accidentales habituales; el control también detecta archivos bloqueados que se
hayan agregado explícitamente a Git.

La configuración de CI, los metadatos exclusivos de GitHub y los archivos locales
excluidos por `.Rbuildignore` no se incorporan al archivo instalable de R.

### Comprobación local

Desde una carpeta de trabajo fuera de la fuente y de los datos privados:

```text
R CMD build /ruta/al/repositorio/varestcf
R CMD check --no-manual varestcf_0.3.0.tar.gz
```

La suite usa inversión densa completa y complemento de Schur como contrastes
independientes. Comprueba alineación, incertidumbre de efectos fijos, covarianzas
cruzadas, varias estructuras residuales y genéticas, proyecciones y contratos de
entrada. Los modelos son sintéticos y pequeños: estas pruebas de exactitud no
certifican capacidad para un número grande de animales. Consulte
[capacidad y límites](CAPACIDAD.md) para planificar recursos.

### Mantenimiento de los dos idiomas

El README y las guías técnicas contienen secciones completas en español e inglés
con enlaces de idioma. La ayuda instalada presenta la referencia en inglés y
una sección en español con argumentos, ecuaciones, salidas y límites. Actualice
ambas versiones cuando cambie un contrato o un modelo admitido. Los nombres de
funciones y argumentos son idénticos en los dos idiomas.

---

<a id="english"></a>

## English

### Repository contents

The repository contains the R reference, help for the two public functions,
documentation, synthetic examples, a consolidated acceptance suite,
authorship/license/citation metadata and GitHub Actions configuration.

Relationships, observations and identifiers for the examples are generated in
`inst/examples/synthetic_models.R`. Dense inversion in that file is an oracle
for small problems, independent of the package's block extraction algorithm.
It is not a recommended method for large models.

Real data, private benchmarks, historical diagnostics and generated files are
kept outside the repository. They are not installed with the package or needed
to run the public tests.

The [first README analysis](../README.md#first-complete-analysis) matches
`inst/examples/getting_started.R`: it can be copied in full or run from the
installed package. Its relationship matrix and observations are explicit; it
does not depend on auxiliary files from the other examples.

### Automatic checks

The [R-CMD-check workflow](../.github/workflows/R-CMD-check.yaml) runs on each
push, pull request and manual dispatch, using R release on Windows and Ubuntu.

1. Reviews paths, formats, sizes and binary content of Git-tracked files using
   the [file guard](../.github/scripts/check-public-files.R).
2. Installs package dependencies and `rcmdcheck`.
3. Builds and installs the package and runs `R CMD check --no-manual`, including
   help, help examples and `tests/acceptance.R`.
4. Saves the results and displays the end of the acceptance suite.

Errors and warnings fail the check. Notes must be reviewed. The numerical suite
runs once per system, inside R CMD check; it is not rerun as an additional step.

External actions are pinned by commit and the workflow uses read permissions.
Its parameters follow the official instructions for
[r-lib/actions](https://github.com/r-lib/actions/tree/v2/check-r-package) and
[actions/checkout](https://github.com/actions/checkout). When updating an action,
review its version and change the pinned commit in the workflow.

These jobs assess portability with the versions installed at each run. They do
not fix a single R/Matrix version for every analysis. Reports record the versions
actually used. This two-system matrix does not certify compatibility with every
R version since the declared minimum.

### Publication review

From the repository root, with the public source tracked or staged in Git:

```text
Rscript .github/scripts/check-public-files.R
```

The guard permits only text formats used by this version, rejects unreviewed
formats and files larger than 1 MiB, and blocks data, outputs, benchmarks and
cache directories. The list is deliberately small: additions requiring another
format must be reviewed before expanding it.

The guard cannot determine whether a matrix written inside a script came from
real data or whether redistribution rights exist. Provenance review is still
required. `.gitignore` prevents common accidental additions; the guard also
detects blocked files explicitly added to Git.

CI configuration, GitHub-only metadata and local files excluded by
`.Rbuildignore` are not included in the installable R archive.

### Local checks

From a working directory outside the source and private data:

```text
R CMD build /ruta/al/repositorio/varestcf
R CMD check --no-manual varestcf_0.3.0.tar.gz
```

The suite uses full dense inversion and the Schur complement as independent
comparisons. It checks alignment, fixed-effect uncertainty, cross-covariances,
multiple residual and genetic structures, projections and input contracts.
Models are small and synthetic: these accuracy tests do not certify capacity
for large animal populations. See [capacity and limits](CAPACIDAD.md#english)
to plan resources.

### Maintaining both languages

The README and technical guides contain complete Spanish and English sections
with language links. Installed function help has an English reference followed
by a Spanish section covering arguments, equations, outputs and limits. Update
both language versions together when changing a contract or supported model.
Function and argument names remain identical in both languages.
