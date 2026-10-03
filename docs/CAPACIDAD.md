# Capacidad de la referencia 0.3

## Alcance disponible

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

## Costos que deben conocerse

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

## Guardas implementadas

`memory_budget_gib` comprueba estimaciones de salida/buffers, inputs, validaciones, ensamblado y almacenamiento tras factorizar C. Las constantes de trabajo incluyen varias copias previstas, pero no garantizan capturar todos los temporales de R/CHOLMOD. Los objetos de entrada ya existen antes de entrar en la función; el presupuesto no controla su creación ni la RAM ocupada por otras aplicaciones.

El ensamblado R puede usar matrices CSC generales, con hasta dos triángulos almacenados temporalmente. Para precisión residual diagonal, se usa una cota estructural a partir de las columnas activas por fila del diseño. Para W general se adopta la cota conservadora M². Una cota excesiva puede rechazar un problema que otro backend admitiría. No debe aumentarse el presupuesto sin evaluar recursos.

La representación actual usa índices CSC de 32 bits: se rechazan cotas de almacenamiento de 2.147.483.647 posiciones o más, y dimensiones fuera del rango admitido. Este límite depende de la representación, no de la RAM disponible. No equivale a un máximo fijo de animales.

El llenado del factor no se conoce con certeza antes de calcularlo. Las guardas son orientativas y portables porque usan dimensiones/objetos, **no una detección garantizada de RAM disponible**. No se promete impedir todas las terminaciones del sistema por falta de memoria. Los residuos normalizados son errores hacia atrás; no garantizan error hacia adelante pequeño si el sistema está mal condicionado.

## Capacidad demostrada y benchmarks pendientes

La suite de aceptación utiliza ocho animales sintéticos con relaciones genealógicas, una corrección genómica inventada, registros repetidos y hasta cuatro coeficientes por animal. Sirve para comprobar exactitud y contratos, no para certificar escalabilidad.

Antes de anunciar capacidades públicas se requiere un benchmark sintético escalonado que registre por separado construcción de Q, validaciones, ensamblado, factorización, extracción y estadístico. Debe variar el número de animales, la densidad de la precisión, los componentes y términos aleatorios, los efectos fijos y los focales, y registrar hardware, hilos, versiones y pico RSS medido externamente. Una única matriz grande no establece una regla universal.

La ruta optimizada futura deberá reproducir esta referencia con tolerancias justificadas. Los límites de recursos no autorizan a modificar las relaciones, omitir efectos presentes en el modelo ni reemplazar covarianzas ausentes por cero.
