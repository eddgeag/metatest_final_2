# METATEST: matemática y relevo para ChatGPT normal

Estado: propuesta metodológica para discusión. No es un método ya aprobado, implementado ni validado en METATEST. No generar código R hasta que el usuario apruebe el diseño.

Este documento separa tres cosas: funcionamiento del framework actual, hallazgos de la auditoría anterior y formulación de una arquitectura candidata. Las cifras de auditoría proceden de los archivos examinados en la conversación anterior; un chat nuevo no debe afirmar que volvió a verificarlas si no recibe y abre esos archivos.

## 1. Objetivo y restricciones

Se buscan firmas biológicas por vista —transcriptómica (TX), proteómica (PR), metabolómica (ME) y clínica (CL)— y una firma predictiva integrada para NP, FD y SO. El predictor utilizará valores observados. Ninguna vista debe tener un cupo obligatorio ni quedar excluida por un ranking que favorezca estructuralmente a las vistas grandes.

| Conjunto |  NP |  FD |  SO | Total |
|----------|----:|----:|----:|------:|
| TRAIN    |   6 |   3 |   7 |    16 |
| TEST     |   4 |   2 |   4 |    10 |

P1, P2 y P3 son análisis científicos independientes. P1 tiene clínica amplia; P2 incluye FLI deliberadamente; P3 excluye FLI y conserva otras clínicas. Sus escenarios MOFA ya se eligieron con TRAIN y deben permanecer congelados:

| Análisis | Escenario                            |
|----------|--------------------------------------|
| P1       | `P1_TODAS_qtx95_qpr50_qme95_qcl40`   |
| P2       | `P2_CON_FLI_qtx95_qpr40_qme95_qcl50` |
| P3       | `P3_SIN_FLI_qtx95_qpr45_qme95_qcl45` |

TEST ya fue inspeccionado durante el desarrollo previo. Desde ahora sirve como auditoría exploratoria externa; no se utiliza para elegir arquitectura, variables, hiperparámetros ni el análisis “ganador”. No adquiere de nuevo condición confirmatoria por cambiar el método.

## 2. Matemática del framework existente

### Reconstrucción MOFA y dependencia entre variables

Para una vista, la parte centrada de la reconstrucción es:

$$
\widehat X^{(v)}=Z(W^{(v)})^T.
$$

Con dos factores, su rango es como máximo dos. Cientos de variables reconstruidas pueden expresar variaciones del mismo espacio de dos dimensiones. Eso no demuestra que sean duplicados exactos, pero impide interpretarlas como cientos de evidencias independientes. Observado y reconstruido provienen de los mismos sujetos: su concordancia tampoco equivale a replicación externa.

### Densidades Student-t descriptivas

Para una variable y una clase, el script estima una localización robusta mediante la mediana y una escala basada en MAD, con alternativas cuando MAD no es utilizable. Regulariza la escala hacia la global:

$$
\omega_k=\frac{n_k}{n_k+\kappa},\qquad
\widetilde\sigma_k^2=\omega_k\sigma_k^2+(1-\omega_k)\sigma_G^2.
$$

La densidad es:

$$
f_k(x)=\frac{1}{\widetilde\sigma_k}
t_4\!\left(\frac{x-\mu_k}{\widetilde\sigma_k}\right),
\qquad
P(Y=k\mid x)=\frac{\pi_k f_k(x)}{\sum_h\pi_h f_h(x)},
\qquad \pi_k=1/3.
$$

Aquí sigma es el parámetro de escala de Student-t, no su desviación estándar. Estas son probabilidades de un modelo univariante ajustado, no probabilidades de que un biomarcador sea verdadero ni una integración bayesiana completa de la incertidumbre de sus parámetros.

La información mutua que calcula el script mediante integración numérica corresponde a:

$$
I(X;Y)=\sum_k\pi_k\int f_k(x)
\log\frac{f_k(x)}{\sum_h\pi_h f_h(x)}\,dx,
\qquad I_{\mathrm{norm}}=\frac{I(X;Y)}{-\sum_k\pi_k\log\pi_k}.
$$

La divergencia Jensen–Shannon de dos clases es:

$$
JS(f_a,f_b)=\tfrac12 KL(f_a\Vert m)+\tfrac12 KL(f_b\Vert m),
\quad m=\tfrac12(f_a+f_b).
$$

El código expresa JS en bits. Las grades combinan estas descripciones con pruebas por permutación, corrección múltiple, patrones y jackknife; no son una escala probabilística calibrada. El framework descriptivo puede conservarse como referencia científica.

### El filtro que se convirtió en universo predictivo

La admisión actual puede escribirse como un producto de indicadores:

$$
G_j=\mathbb 1(A_{rec})\mathbb 1(A_{obs})
\mathbb 1(\text{mismo patrón no incierto})
\mathbb 1(\text{concordancia}\ge2)
\mathbb 1(\text{patrón class-like-high}).
$$

El último conjunto contiene exactamente `NP_like_high`, `FD_like_high` y `SO_like_high`. Una variable con cualquier factor igual a cero no llega al predictor, independientemente de su utilidad conjunta. El problema principal está en convertir una definición de evidencia especialmente fuerte en requisito de elegibilidad predictiva.

## 3. Hallazgos que debe conservar el relevo

### Embudo P3

| Etapa                                    |  TX |  PR |  ME |  CL | Total |
|------------------------------------------|----:|----:|----:|----:|------:|
| Universo                                 | 656 |  13 | 117 |   8 |   794 |
| A reconstruido                           | 631 |  13 | 110 |   7 |   761 |
| Además A observado                       | 229 |   4 |  59 |   1 |   293 |
| Además mismo patrón                      | 221 |   4 |  57 |   1 |   283 |
| Además concordancia ≥2 y class-like-high | 221 |   4 |  57 |   1 |   283 |

La exigencia A observado elimina 468 variables en ese paso, incluidas nueve proteínas y seis clínicas. En P1 el pool final tiene 301 variables; en P2, 290. El cuello de botella inicial no es mRMR.

En P3, permitiendo A/B observado bajo las demás condiciones, el pool sería TX=439, PR=7, ME=96, CL=4: 546 variables. Hasta A/B/C sería 621. Eso ilustra el coste del filtro; no justifica adoptar automáticamente uno de esos universos. T1/T2 tampoco es una solución automática: en P3 reúne 775/794 variables y T2 puede incluir evidencia observada D/E.

### Segundo cuello de botella: ranking y grid

M3-A usa:

$$
R_j=\sum_f |W_{jf}|\,R^2_{f,v(j)}.
$$

M3-B transforma el rango dentro de clase y vista:

$$
R_j^{B}=\frac{\operatorname{rank}(R_j)-0.5}{n_{k,v}}.
$$

Su máximo depende del tamaño de la celda: 0.5 con una variable, 0.875 con cuatro y casi uno en celdas grandes. No garantiza oportunidades comparables entre vistas.

Las grids permiten como máximo siete variables por clase. En P3, las cuatro proteínas A/A tienen rangos por clase A/B de IL-2: 60/22; IL-8: 73/65; VEGF-A: 126/109; Osteonectin: 92/49. PCR queda 56/29. Ninguna llega a ser ajustada. En P2, FLI sobrevive A/A pero queda SO=49/15. Los archivos de exclusión por redundancia no registran descartes en las seis estrategias auditadas.

| Análisis | Firma A: TX/PR/ME/CL | Firma B: TX/PR/ME/CL |
|----------|----------------------|----------------------|
| P1       | 5 / 0 / 11 / 5       | 13 / 0 / 4 / 1       |
| P2       | 8 / 0 / 13 / 0       | 14 / 0 / 4 / 0       |
| P3       | 8 / 0 / 13 / 0       | 14 / 0 / 4 / 0       |

### Otros problemas comprobados

- Cambiar automáticamente a FDR por vista no rescata clínica: en P3, A observado pasa de TX/PR/ME/CL=232/4/63/1 con FDR global a 210/4/88/0 por vista. PCR pasa de q≈0.093 a q≈0.197.
- La CV previa reutiliza pool, ranking y transformaciones aprendidos con todo TRAIN; vuelve a ajustar principalmente el BRMS. No evalúa de forma independiente el procedimiento completo de selección.
- La selección de grid y su evaluación en la misma CV añaden optimismo. La regla llamada 1-SE no prioriza necesariamente la firma más pequeña; P2-A y P3-A seleccionan 21 habiendo opciones menores admisibles. El descarte de configuraciones con clasificación CV perfecta es además una heurística que debe retirarse.
- Con 21 predictores, un multinomial de tres clases estima 42 pendientes y dos interceptos. La convergencia MCMC no demuestra generalización.
- El macro-F1 previo puede omitir una clase fallida al convertir su F1 en NA y promediar con `na.rm=TRUE`. En TEST, P2-B pasa de ≈0.411 a ≈0.274 y P3-B de ≈0.393 a ≈0.262 al incluir correctamente NP con F1=0.
- P2/P3 no son una ablación pura de FLI: cambian también variables clínicas y el universo proteómico. Para aislar su aporte predictivo observado hace falta una comparación pareada con el mismo universo restante, particiones y procedimiento. Eso no mide un efecto causal de FLI.

## 4. Arquitectura candidata: firmas observadas e integración tardía restringida

La propuesta consiste en aprender prototipos regularizados por vista y combinar sus scores mediante un modelo pequeño. No se ha demostrado que supere una alternativa más sencilla. Se plantea como una arquitectura interpretable que merece evaluación, con una referencia multinomial regularizada y un modelo de solo interceptos.

### 4.1 Notación y preprocesamiento

Sea v una de las cuatro vistas, j una variable, k una clase y S el conjunto utilizado para ajustar una etapa concreta. Escribir S importa: dentro de la validación nunca representa automáticamente todo TRAIN.

$$
z_{ivj}=\frac{g_{vj}(x_{ivj})-a_{vj,S}}{h_{vj,S}},\qquad h_{vj,S}>0.
$$

La transformación g, imputación, control de calidad, centrado y escalado se fijan o estiman exclusivamente en S y se aplican sin reajustar al sujeto retenido. Una opción es mediana/MAD, con reglas explícitas para escalas nulas. La elección depende de cómo se procesaron originalmente las plataformas. Las fórmulas siguientes suponen una matriz utilizable tras ese procesamiento.

El universo inicial debe definirse por disponibilidad y calidad técnica dentro del escenario congelado, no por A/A ni por pertenecer a un patrón favorable. Si un filtro supervisado interviene, debe aprenderse dentro del ajuste correspondiente.

### 4.2 Centroides contraídos dentro de cada vista

Con n=\|S\| y K=3, definir:

$$
\mu_{vj,k}=\frac{1}{n_k}\sum_{i\in S:y_i=k}z_{ivj},\qquad
\mu_{vj}=\frac{1}{n}\sum_{i\in S}z_{ivj},
$$

$$
q_{vj}^{2}=\frac{1}{n-K}\sum_k\sum_{i\in S:y_i=k}
(z_{ivj}-\mu_{vj,k})^2,
\qquad
\nu_{vj}=\max\{q_{vj}+\operatorname{mediana}_{\ell}(q_{v\ell}),\epsilon\}.
$$

Aplicar:

$$
m_k=\sqrt{\frac1{n_k}+\frac1n},\quad
d_{vj,k}=\frac{\mu_{vj,k}-\mu_{vj}}{m_k\nu_{vj}},\quad
\widetilde d_{vj,k}=\operatorname{sign}(d_{vj,k})(|d_{vj,k}|-\delta)_+,
$$

$$
\widetilde\mu_{vj,k}=\mu_{vj}+m_k\nu_{vj}\widetilde d_{vj,k}.
$$

Es el componente de centroides contraídos de PAM; delta controla la contracción. Si todos los contrastes de una variable quedan en cero, deja de discriminar. [Tibshirani et al., 2002](https://www.pnas.org/doi/10.1073/pnas.082099299).

La aplicación por vista y la integración posterior son una adaptación propuesta aquí. m_k es el factor de escalado del procedimiento; no debe describirse como el error estándar exacto de una diferencia entre dos medias independientes. La media de clase y la global comparten sujetos.

La firma de una vista es:

$$
\mathcal S_v=\{j:\exists k,\ \widetilde d_{vj,k}\ne0\}.
$$

Los contrastes $\widetilde\mu_{vj,k}-\widetilde\mu_{vj,\ell}$ permiten describir incrementos, disminuciones, patrones compartidos y separadores FD/SO. No se exige una etiqueta única class-like-high. Una firma por vista puede existir aunque esa vista no añada información al integrador.

Para empezar se plantea un delta común, con una decisión explícita sobre su valor o un conjunto pequeño de candidatos. No se propone buscar libremente cuatro umbrales independientes. Con FD muy pequeño, el mismo desplazamiento de media recibe más contracción; protege frente a estimaciones inestables, pero también puede ocultar señales FD reales. La escala común de d no neutraliza por sí sola la multiplicidad de una vista grande.

### 4.3 De la firma a scores de clase

Para un sujeto nuevo, definir por vista:

$$
r_{vk}(z)=-\frac12\sum_j
\frac{(z_{vj}-\widetilde\mu_{vj,k})^2}{\nu_{vj}^2},
\qquad
u_{vk}(z)=r_{vk}(z)-\frac1K\sum_h r_{vh}(z).
$$

Una normalización candidata, calculada en el mismo S, es:

$$
c_{v,S}=\max\left\{\epsilon,
\sqrt{\frac1{nK}\sum_{i\in S}\sum_k u_{vk}(z_i)^2}\right\},
\qquad s_{vk}(z)=\frac{u_{vk}(z)}{c_{v,S}}.
$$

Hay tres scores centrados por vista, con solo dos contrastes independientes. No se añaden priors de clase a cada vista para volver a contarlos cuatro veces; el componente de clase se introduce en el integrador.

Esta normalización cambia la escala, pero no elimina la redundancia ni demuestra equidad entre tamaños de vista. También puede amplificar un bloque de señal débil. Si todos los contrastes son cero, los scores son cero. Antes de implementar debe resolverse si esta escala, una alternativa más conservadora o una escala prefijada resulta defendible; no explorar muchas versiones y escoger la que mejor se vea.

La métrica diagonal no estima las correlaciones entre variables. Duplicados y medidas casi equivalentes requieren reglas de calidad o sensibilidad explícitas. Aprender redes grandes de correlación con N=16 añadiría otra fuente de inestabilidad.

### 4.4 Integrador multinomial restringido

$$
\eta_{ik}=\alpha_k+\sum_{v=1}^{4}b_v s_{ivk},\qquad
p_{ik}=\frac{\exp(\eta_{ik})}{\sum_h\exp(\eta_{ih})},
\qquad \alpha_{NP}=0,\quad b_v\ge0.
$$

Cada vista tiene un único peso compartido por las clases. Hay cuatro pesos y dos interceptos en esta etapa. Eso no significa seis parámetros aprendidos en todo el procedimiento: también se estiman preprocesamiento, centroides, escalas y selección.

No imponer $\sum_v b_v=1$: debe ser posible que todas las vistas reciban peso cero y sobreviva el modelo de interceptos. Compartir b_v reduce flexibilidad, pero puede desaprovechar una vista útil para un contraste e inútil para otro. Los pesos tampoco son efectos causales ni medidas estables de importancia cuando las vistas son redundantes.

Una función objetivo posible es:

$$
\min_{\alpha,\,b\ge0}
-\frac1n\sum_{i\in S}\log p_{i,y_i}
+\lambda\left[\rho\sum_v b_v+
\frac{1-\rho}{2}\sum_v b_v^2\right]
+\frac{\lambda_{\alpha}}2\sum_{k\ne NP}(\alpha_k-\alpha_{0k})^2.
$$

Con 0\<rho\<1 prefijado, la parte L1 permite pesos exactamente cero y la parte cuadrática estabiliza. lambda_alpha y el centro alpha_0 deben declararse. Un centro cero regulariza hacia igualdad; otro basado en frecuencias exige justificar el objetivo. Las proporciones de esta muestra no constituyen prevalencias poblacionales demostradas. Minimizar esta función no equivale a obtener una distribución posterior bayesiana completa.

### 4.5 Entrenamiento del integrador sin scores autoajustados

El integrador debe entrenarse con scores obtenidos dejando fuera al sujeto o grupo que recibe el score:

$$
s^{CF}_{ivk}=s^{\,S\setminus F(i)}_{vk}(x_i).
$$

F(i) es el fold interno que contiene i. Se recalculan preprocesamiento, selección, centroides y escala en su complemento. Estos scores cruzados se utilizan para ajustar alpha y b.

La pérdida calculada sobre esos mismos sujetos después de ajustar el integrador sigue siendo una pérdida de entrenamiento del integrador. No es la evaluación externa de toda la arquitectura.

### 4.6 Contribución exacta de las variables

Para un modelo ya ajustado, la diferencia de distancias cuadráticas cancela los términos z². Por tanto:

$$
\log\frac{p_k(z)}{p_{NP}(z)}
=A_k+\sum_{v,j}\beta^{eff}_{vj,k}z_{vj},
$$

$$
\boxed{\beta^{eff}_{vj,k}=
\frac{b_v}{c_{v,S}\nu_{vj}^{2}}
\left(\widetilde\mu_{vj,k}-\widetilde\mu_{vj,NP}\right)}.
$$

El intercepto incluye:

$$
A_k=\alpha_k-\alpha_{NP}
-\frac12\sum_{v,j}\frac{b_v}{c_{v,S}\nu_{vj}^2}
\left(\widetilde\mu_{vj,k}^{\,2}-\widetilde\mu_{vj,NP}^{\,2}\right).
$$

Esto permite informar la contribución aditiva $\beta^{eff}_{vj,k}z_{ivj}$ a los log-odds de cada sujeto. Si g es identidad, el coeficiente respecto a x es beta_eff/h; si g es logarítmica, la linealidad es respecto a la variable transformada. La interpretación depende de la referencia y del centrado. No mide causalidad y puede repartirse entre variables correlacionadas.

## 5. Cómo incorporar evidencia sin otro filtro binario

La versión inicial más auditable mantiene las grades congeladas como anotación. Se informa qué miembros de cada firma tienen A/A, B, C, concordancia o señal reconstruida, sin permitir que esas etiquetas calculadas con todo TRAIN decidan la elegibilidad dentro de la CV.

Una extensión sería:

$$
\widetilde d_{vj,k}=\operatorname{sign}(d_{vj,k})
(|d_{vj,k}|-\delta e_{vj})_+,\qquad e_{vj}>0.
$$

Una evidencia más débil implicaría mayor e y mayor contracción. Pero fijar A=1, B=2 y C=4 sería una decisión adicional, no una calibración estadística conocida. Si e se aprende de estos sujetos, debe recalcularse dentro de cada ajuste. Si depende de reconstrucción MOFA, también debe especificarse qué se reajusta y qué permanece condicionado al escenario congelado. Con un solo FD en ajustes internos, ciertas estadísticas pueden no ser estimables; no se sustituyen silenciosamente por las calculadas con todo TRAIN.

Esta extensión queda pendiente. El modelo básico usa e=1 y separa evidencia descriptiva de aprendizaje predictivo.

## 6. Validación con tres sujetos FD

Un diseño candidato es CV externa estratificada de tres folds: cada validación contiene un FD y cada entrenamiento externo, dos. Dentro del entrenamiento externo puede utilizarse leave-one-out para construir los scores cruzados del integrador. Al retirar un FD, el modelo base tendrá solo uno. La varianza agrupada puede calcularse si hay grados de libertad en las demás clases, pero el centro FD será extremadamente incierto.

No usar “3×3 anidada estratificada” como receta automática: con dos FD en el entrenamiento externo no se puede repartir un FD en cada uno de tres folds internos.

Para cada fold externo:

1.  Apartar sus sujetos antes de aprender cualquier componente del nuevo predictor.
2.  Construir scores cruzados dentro del entrenamiento externo y ajustar el integrador con ellos.
3.  Reajustar las firmas de base usando todo ese entrenamiento externo.
4.  Aplicar bases e integrador al fold externo, sin reajustes sobre sus sujetos.

Este procedimiento presupone hiperparámetros fijados. Si delta o lambda se eligen mediante validación, la selección debe hacerse exclusivamente dentro del entrenamiento externo y evaluar allí el procedimiento de dos etapas completo. No seleccionar lambda por la pérdida de entrenamiento del integrador y presentarla como CV. La anidación adicional puede dejar clases sin representación; con estos datos, esa limitación puede obligar a prefijar parámetros o preferir un método de una sola etapa.

La selección de modelos puede sobreajustar su propio criterio de evaluación. [Cawley y Talbot, 2010](https://jmlr.org/papers/v11/cawley10a.html).

Dos límites permanecen aunque se corrija la validación del predictor:

- Los escenarios ya se seleccionaron con todo TRAIN. Mantenerlos congelados implica evaluar condicionalmente a esa selección previa; no permite afirmar una estimación independiente de todo el proceso de descubrimiento.
- Repetir particiones reutiliza los mismos 16 sujetos. Su dispersión describe sensibilidad a las particiones; dividirla por la raíz del número de repeticiones no proporciona automáticamente un intervalo de incertidumbre poblacional. [Varoquaux, 2018](https://arxiv.org/abs/1706.07581).

Después de fijar el diseño, ajustar las bases con todo TRAIN y el integrador con scores cruzados de TRAIN. Aplicar una vez a TEST para la auditoría prevista, sin retocar el método por sus resultados. No garantizar calibración de probabilidades con esta muestra.

## 7. Métricas, estabilidad y aporte de una vista

Usar predicciones externas al ajuste de ambas etapas. Agregar los sujetos retenidos antes de calcular métricas de una repetición:

$$
LL=-\frac1n\sum_i\log p^{out}_{i,y_i},\qquad
Brier=\frac1n\sum_i\sum_k(p^{out}_{ik}-\mathbb1[y_i=k])^2,
$$

$$
BA=\frac1K\sum_k\frac{TP_k}{TP_k+FN_k},\qquad
F1_{macro}=\frac1K\sum_k\frac{2TP_k}{2TP_k+FP_k+FN_k}.
$$

El Brier aquí no se divide por K y tiene rango 0–2. Si una clase tiene soporte real y TP=0, su F1 es cero. Informar pérdida y sensibilidad por clase; macro-AUC puede describirse con cautela, pero no debe dominar el tuning con FD=3. El logloss medio da más peso a las clases frecuentes: declarar si el objetivo principal es la distribución observada o una ponderación igual por clase. No cambiarlo tras ver resultados.

Para evaluar aporte incremental de una vista, reajustar el procedimiento sin ella sobre las mismas particiones:

$$
\Delta_v=LL_{-v}-LL_{todas}.
$$

Un valor positivo favorece su contribución predictiva incremental bajo este procedimiento. Un valor próximo a cero no demuestra ausencia de biología: puede reflejar redundancia o escasa potencia. No buscar libremente entre todas las combinaciones de vistas.

Para estabilidad descriptiva:

$$
\widehat f_{vj}=\frac1R\sum_{r=1}^{R}
\mathbb1(j\in\mathcal S_v^{(r)}).
$$

Registrar también consistencia de signo, contraste, peso de vista y sensibilidad al retirar cada FD. Las frecuencias de selección no son probabilidades posteriores ni control formal de falsos descubrimientos. Si se convierten en un nuevo filtro para predecir, ese filtro también debe aprenderse dentro de la validación.

## 8. Alternativas que deben compararse críticamente

| Alternativa | Formulación o principio | Límite relevante en METATEST |
|----|----|----|
| Multinomial ridge/elastic net | Minimizar logloss más penalización de pendientes observadas | Una sola etapa facilita validación; el tamaño y correlación de las vistas siguen importando |
| Sparse-group multinomial | Penalización de grupos de coeficientes por vista y de elementos individuales | Los pesos de grupo no deben elegirse arbitrariamente; una corrección por tamaño no garantiza igualdad de oportunidades |
| Bayesiano jerárquico | Escalas de contracción global, por vista y por variable | Muchos parámetros débilmente identificados; los priors pueden dominar con N=16 |
| DIABLO/block-sPLS-DA | Componentes supervisados que combinan discriminación y relaciones entre bloques | Número de componentes, sparsidad y diseño entre bloques amplían el tuning |
| Pareto | Conservar soluciones no dominadas en predicción, tamaño y estabilidad | Organiza compromisos; no es por sí mismo un predictor ni evita sobreajuste de selección |

Una expresión genérica de sparse-group multinomial, con B_v la matriz de pendientes de una vista, es:

$$
\min_{\alpha,B}\mathcal L(\alpha,B)
+\lambda\left[(1-\rho)\sum_v w_v\lVert B_v\rVert_F
+\rho\sum_{v,j,k}|B_{vj,k}|\right].
$$

Una jerarquía bayesiana esquemática sería:

$$
\beta_{vj,k}\mid\tau,\tau_v,\lambda_{vj}
\sim N(0,\tau^2\tau_v^2\lambda_{vj}^2).
$$

Esta última ecuación está incompleta sin priors, identificación y regularización de efectos grandes: no debe etiquetarse automáticamente como regularized horseshoe. Los priors continuos tampoco producen ceros exactos como L1. Estas alternativas necesitan su propia especificación; no deben sumarse todas a una búsqueda amplia con 16 sujetos.

La integración por centroides se propone por trazabilidad y reducción de parámetros en la última etapa. Su superioridad, estabilidad, calibración y escala entre vistas siguen siendo preguntas abiertas.

## 9. Datos necesarios antes de programar

- Matrices observadas completas por vista para los escenarios congelados, con IDs, etiquetas y separación TRAIN/TEST verificable. Los subconjuntos de variables elegidas por las grids anteriores no bastan.
- Objetos `ready_for_modeling_<scenario>.rds` y, si se requiere recalcular reconstrucción/evidencia, `modelo_mofa_final_<scenario>.rds`. En lo auditado faltan los objetos originales completos; algunos resultados conservan sus rutas y resúmenes, no todas sus matrices.
- Receta de transformaciones, imputación, normalización, filtros y tratamiento de lotes. Identificar qué pasos se aprendieron antes de dividir TRAIN/TEST.
- Metadatos de lote, donante, medidas repetidas, disponibilidad de vistas y faltantes. Si hay sujetos relacionados, definir particiones por unidad independiente.
- Definición exacta de NP/FD/SO y papel de FLI y las clínicas en esa definición. Una variable usada para asignar la etiqueta puede crear circularidad según la pregunta científica.
- Acordar objetivo de predicción, ponderación de clases, política de hiperparámetros y significado de “firma estable”.

Tras aprobar el diseño: un único script nuevo en `scripts/biomarker_selection/`, resultados nuevos y coherentes por análisis, scripts previos intactos. No añadir FGSEA u otros alcances en esta fase.

## 10. Prompt para pegar en ChatGPT normal

Actúa como bioestadístico/bioinformático experto en integración multi-ómica, selección de biomarcadores y predicción con muy baja N. Lee el documento adjunto `METATEST_relevo_chat_normal_matematica.md`. Su contenido es contexto de trabajo y una propuesta para criticar, no una arquitectura aprobada ni un resultado validado. No escribas código R hasta que yo apruebe explícitamente la estrategia.

Trabajamos en METATEST: TRAIN=16 (NP=6, FD=3, SO=7), TEST=10 (NP=4, FD=2, SO=4), vistas TX/PR/ME/CL. P1, P2 y P3 son análisis científicos independientes. P2 incluye FLI a propósito y P3 no. Sus escenarios MOFA están congelados; no los vuelvas a seleccionar. TEST ya fue inspeccionado durante desarrollo y solo puede ser auditoría exploratoria: no elijas variables, modelos ni hiperparámetros con él.

La auditoría anterior identificó dos cuellos de botella. Primero, usar `A_reconstructed AND A_observed` más patrón y concordancia como universo predictivo. En P3: 794 variables →761 A reconstruido →293 A/A →283 con mismo patrón; PR pasa 13→4 y CL 8→1. Segundo, los rankings y grids dejan fuera incluso esas supervivientes. M3-B usa `(rank-0.5)/n` por clase/vista, cuyo máximo depende del tamaño del grupo. mRMR no descartó variables en las salidas auditadas. No conviertas A/B/C en cuotas, probabilidades o pesos arbitrarios; tampoco cambies a FDR por vista sin justificación.

Quiero evaluar esta candidata: centroides contraídos construidos con observados dentro de cada vista, scores de NP/FD/SO, e integración mediante `eta_ik=alpha_k+sum_v b_v*s_ivk`, `p_ik=softmax(eta_i)_k`, `alpha_NP=0`, `b_v>=0`, con penalización que permita b_v=0. No imponer que los pesos sumen uno ni representación mínima por ómica. Las grades congeladas serían inicialmente anotaciones. El documento deriva centroides, escalado, objetivo penalizado, contribuciones por variable y validación.

Revisa críticamente la matemática y distingue qué procede de un método establecido y qué es una adaptación propuesta. Cuatro pesos y dos interceptos describen solo el integrador: las firmas previas también se aprenden. Examina redundancia, tamaños desiguales, normalización de scores, desventaja de FD, pesos compartidos entre clases y el riesgo de amplificar ruido. Si una alternativa más sencilla resulta más defendible, arguméntalo.

Toda selección, transformación y construcción de scores debe ocurrir dentro del ajuste correspondiente. Los scores cruzados usados para entrenar el integrador no son por sí solos predicciones externas del procedimiento entero. Explica una validación viable con FD=3; no propongas una anidación que deje clases vacías. No presentes repetición de particiones como aumento de N. La selección previa del escenario con todo TRAIN limita lo que puede afirmarse de la validación posterior.

Compara la candidata con una referencia multinomial ridge/elastic net y al menos otra opción razonable, sin diseñar una búsqueda masiva de modelos. Separa evidencia biológica, utilidad predictiva marginal y aporte incremental al combinar vistas. Corrige el macro-F1 para incluir clases con cero aciertos. P2/P3 no aíslan únicamente FLI: cualquier ablación debe usar el mismo universo restante y particiones, sin interpretación causal.

En tu primera respuesta: explica con ecuaciones qué conservarías o cambiarías; identifica las decisiones todavía abiertas; recomienda una arquitectura provisional con sus límites; concreta cómo validarla; y enumera únicamente los datos faltantes necesarios para implementarla. No repitas toda la auditoría ni digas que inspeccionaste archivos que no he adjuntado. Pregunta por información solo cuando cambie una decisión relevante. Después discutiremos la estrategia y, únicamente tras mi aprobación, producirás un único script nuevo sin modificar los anteriores.

Si adjunto los dos scripts R y `lo.zip`, úsalos para verificar los puntos que lo requieran. Si no están adjuntos, trata las cifras de este documento como hallazgos reportados de la auditoría anterior, no como una nueva verificación tuya.
