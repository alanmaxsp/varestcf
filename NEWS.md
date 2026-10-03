# 0.3.0 — varios términos aleatorios

- Primera publicación pública el 2026-10-03: instalación fija por `v0.3.0`, documentación bilingüe y 147 comprobaciones sintéticas aprobadas en Windows y Linux (first public release with pinned installation, bilingual documentation and 147 passing synthetic checks on Windows and Linux).
- Validación explícita de rango fijo mediante QR dispersa escalada para rechazar diseños no identificados de forma consistente entre sistemas (explicit fixed-design rank validation across backends).
- Especificación de términos por nombre, identificación, base y covarianza.
- Selección explícita del efecto genético focal; todos los términos permanecen en las ecuaciones.
- Covarianzas conjuntas entre términos que comparten una matriz de relaciones.
- Repetibilidad, efectos maternos, efectos comunes y regresiones aleatorias con términos adicionales.
- Diseño fijo y precisión residual compartidos durante la preparación.
- Controles de equivalencia con ecuaciones explícitas y contratos de entrada.
- Alcance en modelos mixtos lineales gaussianos.

# 0.2.1 — alcance general

- Documentación y ejemplos organizados por estructuras matemáticas del modelo.
- Cobertura residual, términos aleatorios adicionales y prioridades de ampliación explícitas.
- Controles sintéticos de repetibilidad, efectos maternos y heterogeneidad por clases.

# 0.2.0 — interfaz automática local

- Datos y fórmulas, con A/H o su inversa suministrada; sin construcción de relaciones.
- Diseños, precisión conjunta, alineación y EBV compatibles preparados automáticamente.
- Regresión aleatoria con proyecciones completas y múltiples respuestas con ausentes.
- Se conserva la interfaz matricial y las dos funciones públicas.

# 0.1.0 — referencia local

- Dos funciones exportadas: reconstrucción focal completa y VarEst-CF.
- Diseños y precisión conjunta explícitos, con alineación por identificadores.
- Validación de matrices, incertidumbre fija, residuos y covarianzas focales.
- Diagnósticos de recursos y tiempos por etapa; límites conservadores documentados.
- Ejemplos sintéticos de efectos animales, múltiples caracteres y regresión aleatoria; una suite consolidada de aceptación.
- Sin dependencia de BLUPF90/Rcpp, datos reales, caches ni modalidad de PEV/PEC nativa.
