# Cobertura de modelos y prioridades

La unidad de diseño es un modelo mixto lineal y un objetivo genético focal. El nombre de una aplicación no determina la implementación.

## Qué se puede usar actualmente

| Estructura | Interfaz de fórmulas | Interfaz matricial |
|---|---|---|
| Un carácter y un efecto genético con relaciones suministradas | Disponible; K o K⁻¹ y varianza genética | Disponible |
| Varios caracteres, covarianzas genéticas y observaciones ausentes | Disponible; mismos términos fijos, coeficientes separados por carácter | Diseños explícitos pueden diferir entre caracteres |
| Regresión aleatoria para un carácter | Disponible; base especificada con fórmula y covarianza de coeficientes | Disponible |
| Residuos heterogéneos para una respuesta, incluidos modelos de regresión aleatoria | Vector por fila o expresión conocida; admite clases mediante un vector | Precisión diagonal por observación |
| Residuos correlacionados entre caracteres de una misma fila | Matriz residual común; se invierte el subbloque observado cuando hay ausentes | Precisión residual completa |
| Residuos heterogéneos por clase en modelos multicaracter o correlaciones entre filas | No hay entrada automática para estas estructuras | Matriz de precisión residual completa |
| Repetibilidad con ambiente permanente | Disponible con términos separados; contrastado con ecuaciones explícitas | Disponible |
| Efectos genéticos directos y maternos correlacionados | Disponible con covarianza conjunta y una matriz de relaciones común; contrastado | Disponible |
| Camada, ambiente común o grupos aleatorios adicionales | Disponible; identificación por niveles y varianza/covarianza suministrada; ejemplo de grupo contrastado | Disponible |
| Regresión aleatoria para una respuesta, con regresión de ambiente permanente | Disponible; bases y covarianzas por término, residuos heterogéneos; contrastado | Disponible |
| Términos de regresión aleatoria correlacionados entre sí | Disponible con covarianza conjunta y relaciones comunes; ejemplo contrastado | Disponible |
| Varios caracteres con términos aleatorios adicionales | Disponible con intercepto por carácter; caso con ambiente permanente, residuos correlacionados y ausentes contrastado | Disponible |
| Regresión aleatoria multicaracter | Pendiente de preparación automática | Representable con diseños y precisión conjunta; sin ejemplo dedicado actual |
| Dominancia, efectos sociales u otros kernels | Requieren términos y objetivo genético explícitos | Representables si la covarianza conjunta es definida positiva; sin ejemplos dedicados actuales |
| Precisión singular, grupos genéticos o parametrizaciones con restricciones especiales | Sin soporte específico | Debe suministrarse una parametrización identificada con precisión definida positiva |

“Representable” significa compatible con la formulación lineal general, no una garantía de validación de cada familia ni de capacidad computacional para cualquier tamaño. No se estiman componentes de varianza.

El alcance se limita a modelos mixtos lineales gaussianos.

## Por qué la estructura residual importa

Sea y = Xb + Zu + e, con Cov(u) = B, Cov(e) = R y u independiente de e. Defina Q = B⁻¹ y W = R⁻¹. Las ecuaciones completas son

```text
C = [ X'WX       X'WZ       ]
    [ Z'WX       Z'WZ + Q  ]
C [b_hat; u_hat] = [X'Wy; Z'Wy]
Sigma_f = E_f' C⁻¹ E_f
```

Con residuos independientes y heterogéneos, R = diag(v₁,…,vₙ) y W = diag(1/v₁,…,1/vₙ). Las varianzas deben ser positivas y estar asociadas a las filas originales correctas. Cambiar vᵢ cambia C, las predicciones y Sigma_f; no basta con corregir el estadístico al final. Varianza heterogénea y correlación residual son estructuras distintas.

Para regresión aleatoria, u_i(t) = phi(t)'a_i. Si Sigma_ab es el bloque de errores entre los coeficientes a y b de los focales,

```text
Sigma_f(t,s) = sum_ab phi_a(t) phi_b(s) Sigma_ab
u_hat_f(t)   = sum_a  phi_a(t) a_hat_f,a
```

La evaluación actual calcula Sigma_f(t,t) para cada punto pedido y aplica VarEst-CF a los animales focales en ese punto. No devuelve automáticamente una matriz entre diferentes puntos.

Con P = I − 11'/n y denominador d = n o n−1:

```text
VarEst-CF = [u_hat_f' P u_hat_f + trace(P Sigma_f)] / d
```

R determina la incertidumbre de predicción; el objetivo sigue siendo la variación del componente genético focal, no la varianza residual ni la varianza fenotípica total.

## Efectos que deben permanecer en las ecuaciones

En un modelo de repetibilidad, y = Xb + Z_a a + Z_p p + e. Aunque sólo interese a, p debe incluirse con su covarianza. Use Z = [Z_a Z_p] y la precisión conjunta de (a,p). Seleccione como focales únicamente los coeficientes de a. El bloque focal se extrae del sistema completo; no se obtiene eliminando p y reconstruyendo un modelo reducido.

El mismo criterio vale para efectos maternos y cualquier otro término. Las covarianzas entre efectos directos y maternos deben quedar en la covarianza conjunta antes de invertirla. Para objetivos combinados se necesita además la transformación lineal de las predicciones y de todos los bloques de error correspondientes. La interfaz actual no automatiza objetivos genéticos arbitrarios.

## Especificación de varios términos

Use `random_effects` como lista nombrada. Cada término tiene una columna `id`, una base `random` (predeterminada: intercepto) y una `variance`. Puede suministrarse una matriz de relaciones o su inversa; si se omite, los niveles son independientes. `levels` permite incluir niveles independientes no observados; en matrices suministradas se conservan todas las filas/columnas.

Seleccione un único término genético mediante `focal_effect`. Sólo sus coeficientes focales se extraen, conservando todos los demás términos en el sistema. Los efectos no se identifican como genéticos por el nombre: el usuario define el objetivo y proporciona su matriz de relaciones.

Los grupos de `correlated_effects` incluyen una lista ordenada de nombres y una covarianza conjunta. Sus términos deben compartir la misma matriz de relaciones. No se admiten grupos superpuestos. Si se declaran varianzas individuales además de la conjunta, sus bloques diagonales deben coincidir.

La covarianza conjunta sigue el orden de los términos y luego el de sus componentes. En una matriz nombrada, un término escalar usa su nombre; un término con varios componentes usa `término:componente`, por ejemplo `direct:(Intercept)` y `direct:time`. La función alinea los nombres y las IDs de relaciones antes de construir los bloques.

Las IDs faltantes en filas observadas causan un error. `allow_missing = TRUE` es una declaración explícita de incidencia cero de ese término. No debe usarse para omitir un efecto desconocido que el modelo requiere: represente ese nivel y su covarianza.

## Ampliaciones pendientes dentro del alcance lineal

1. Regresión aleatoria multicaracter.
2. Estructuras residuales multicaracter por clase y correlaciones entre registros; diseños fijos específicos por carácter.
3. Transformaciones lineales generales del objetivo focal, incluida combinación de componentes.

La ampliación debe conservar dos funciones públicas y una referencia matricial transparente. Cada término necesita una columna de identificación, una base de diseño y una covarianza suministrada; el usuario debe distinguir los términos presentes en el modelo del objetivo que quiere resumir.

## Fuentes sobre estructuras relevantes

La [documentación de modelos animales de Genstat](https://genstat.kb.vsni.co.uk/knowledge-base/animal-model/) describe efectos maternos, ambientales comunes y términos aleatorios adicionales. El [tutorial de Yutaka Masuda sobre heterogeneidad residual](https://masuday.github.io/blupf90_tutorial/vc_advanced_aireml.html) presenta varianzas por clases y funciones de covariables. Su [índice de modelos de evaluación genética](https://masuday.github.io/blupf90_tutorial/index.html) incluye repetibilidad, modelos maternos, regresión aleatoria y dominancia. Estas fuentes fundamentan las estructuras a considerar; no certifican el soporte de este paquete.
