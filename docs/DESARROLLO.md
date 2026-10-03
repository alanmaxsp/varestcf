# Desarrollo y comprobaciones

## Contenido del repositorio

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

## Chequeo automático

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

## Revisión de publicación

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

## Comprobación local

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
