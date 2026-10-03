# Desarrollo y comprobaciones / Development and checks

<p align="center">
  <a href="#espanol">Español</a> · <a href="#english">English</a>
</p>

<a id="espanol"></a>

## Español

### Ramas y publicación de versiones

| Referencia | Uso |
|---|---|
| `main` | Línea estable para uso habitual; sólo recibe cambios revisados |
| `develop` | Cambios y pruebas para la próxima versión |
| `vX.Y.Z` | Etiqueta de una versión publicada que se conserva sin cambios |

El nombre de una rama no certifica por sí mismo la calidad del código. Esta
separación es una política de trabajo; no se han configurado reglas de protección
que impidan escribir directamente en `main`. Mantenga `main` como rama
predeterminada del repositorio.

Para el trabajo habitual en GitHub Desktop:

1. Abra o clone el repositorio y use **Fetch origin** para consultar novedades.
   Aplique **Pull origin** si hay cambios pendientes. Antes de cambiar de rama,
   guarde sus cambios mediante un commit en su rama de trabajo.
2. En **Current branch**, seleccione `develop`. Para una tarea aislada puede crear
   una rama nueva desde `develop`, por ejemplo `feature/nueva-funcion`, e
   incorporarla posteriormente a `develop` mediante una pull request.
3. Edite y pruebe allí. Revise **Changes**, haga el commit y use **Push origin**.
   Esto actualiza desarrollo; no actualiza `main` ni versiones ya publicadas.
4. Cuando el cambio esté listo, abra una **pull request** con destino `main` y
   origen `develop`. Revise el diff y los resultados de las comprobaciones.
   Las pruebas aprobadas no sustituyen la revisión metodológica.
5. Incorpore la pull request sólo después de revisar esos resultados. Esto cambia
   lo que reciben quienes instalan desde `main`. **No elimine `develop`** después
   de incorporarla: es la rama de desarrollo permanente.
6. Para publicar una versión, actualice `DESCRIPTION`, `NEWS.md` y las referencias
   de instalación/cita necesarias dentro del cambio revisado. En **Releases**,
   cree una publicación con etiqueta `vX.Y.Z` sobre el commit validado de `main`.
   No mueva ni reemplace etiquetas anteriores; publique una nueva para corregirlas.
7. Después de publicar, sincronice `develop` incorporando los cambios de `main`
   sin sobrescribir trabajos pendientes. Durante desarrollo puede usar una versión
   como `X.Y.Z.9000` en `DESCRIPTION`; la versión publicada debe volver a una
   numeración de publicación apropiada.

Una release identifica y describe una versión del código; no implica publicación
en CRAN ni creación automática de un instalador binario de R. Crear `develop`
tampoco crea una release. Antes de una primera publicación numerada, las
instalaciones por etiqueta que muestra el README son ejemplos pendientes.

Las ramas y etiquetas pertenecen al mismo repositorio. Cuando éste es público,
desarrollo también lo es. La separación controla qué versión se instala por
defecto, no oculta el código en desarrollo. Instalar otra rama en la misma
biblioteca de R reemplaza el paquete instalado; para comparar instalaciones
simultáneas use bibliotecas separadas.

Fuentes: [ramas de GitHub](https://docs.github.com/en/pull-requests/reference/branches),
[releases de GitHub](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases)
y [selección de referencias en remotes](https://remotes.r-lib.org/reference/install_github.html).

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

### Branches and version releases

| Reference | Purpose |
|---|---|
| `main` | Stable line for normal use; receives reviewed changes only |
| `develop` | Changes and tests for the next version |
| `vX.Y.Z` | Published version tag retained unchanged |

A branch name does not certify code quality by itself. This separation is a
working policy; no protection rules have been configured to prevent direct writes
to `main`. Keep `main` as the repository's default branch.

For routine work in GitHub Desktop:

1. Open or clone the repository and use **Fetch origin** to check for updates.
   Use **Pull origin** if updates are pending. Before switching branches, save
   your changes in a commit on your working branch.
2. Under **Current branch**, select `develop`. For an isolated task, create a new
   branch from `develop`, such as `feature/new-function`, and later merge it into
   `develop` through a pull request.
3. Edit and test there. Review **Changes**, commit and use **Push origin**. This
   updates development; it does not update `main` or published versions.
4. When ready, open a **pull request** with `main` as base and `develop` as head.
   Review the diff and check results. Passing tests do not replace methodological
   review.
5. Merge only after reviewing those results. This changes what users installing
   from `main` receive. **Do not delete `develop`** after merging: it is the
   permanent development branch.
6. To release a version, update `DESCRIPTION`, `NEWS.md` and any required
   installation/citation references within the reviewed change. Under
   **Releases**, create a release with tag `vX.Y.Z` on the validated `main` commit.
   Do not move or replace previous tags; publish a new version to fix them.
7. After release, synchronize `develop` by merging changes from `main` without
   overwriting pending work. During development, a version such as `X.Y.Z.9000`
   may be used in `DESCRIPTION`; the published version must return to an
   appropriate release version number.

A release identifies and describes a code version; it does not imply CRAN
publication or automatically create an R binary installer. Creating `develop`
does not create a release either. Until the first numbered release, the README's
tag installation commands are pending examples.

Branches and tags belong to the same repository. When it is public, development
is public too. This separation controls the default installation version; it
does not hide development code. Installing another branch in the same R library
replaces the installed package; use separate libraries to compare simultaneous
installations.

Sources: [GitHub branches](https://docs.github.com/en/pull-requests/reference/branches),
[GitHub releases](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases)
and [reference selection in remotes](https://remotes.r-lib.org/reference/install_github.html).

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
