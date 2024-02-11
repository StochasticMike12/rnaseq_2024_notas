################################################################################
#                                                                              #
#     5. Modelos estadísticos                                 01/02/2024       #
#                                                                              #
################################################################################

## ?model.matrix
mat <- with(trees, model.matrix(log(Volume) ~ log(Height) + log(Girth)))
mat

colnames(mat)

summary(lm(log(Volume) ~ log(Height) + log(Girth), data = trees))

## Datos de ejemplo
(sampleData <- data.frame(
  genotype = rep(c("A", "B"), each = 4),
  treatment = rep(c("ctrl", "trt"), 4)
))

## Creemos las imágenes usando ExploreModelMatrix
vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment,
  textSizeFitted = 4
)

## Veamos las imágenes
cowplot::plot_grid(plotlist = vd$plotlist)

## Cargar ExploreModelMatrix para poder usar la función ExploreModelMatrix
library(ExploreModelMatrix)

## Usaremos shiny otra ves
app <- ExploreModelMatrix(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment
)
if (interactive()) shiny::runApp(app)

## Valor esperado = media = promedio
## Indicator = 0 si es el valor de referencia
## Indicator = 1 si no es el valor de referencia
## E[y|ctr~gen]-E[y|trt~gen]

## EJEMPLO 2
(sampleData <- data.frame(
  Response = rep(c("Resistant", "Sensitive"), c(12, 18)),
  Patient = factor(rep(c(1:6, 8, 11:18), each = 2)),
  Treatment = factor(rep(c("pre","post"), 15)),
  ind.n = factor(rep(c(1:6, 2, 5:12), each = 2))))

vd <- VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ Response + Response:ind.n + Response:Treatment,
  textSizeFitted = 3
)
cowplot::plot_grid(plotlist = vd$plotlist, ncol = 1)

app <- ExploreModelMatrix(
  sampleData = sampleData,
  designFormula = ~ Response + Response:ind.n + Response:Treatment
)

if (interactive()) shiny::runApp(app)


## EJEMPLO 3

(sampleData = data.frame(
  condition = factor(rep(c("ctrl_minus", "ctrl_plus",
                           "ko_minus", "ko_plus"), 3)),
  batch = factor(rep(1:6, each = 2))))

vd <- VisualizeDesign(sampleData = sampleData,
                      designFormula = ~ 0 + batch + condition,
                      textSizeFitted = 4, lineWidthFitted = 20,
                      dropCols = "conditionko_minus")
cowplot::plot_grid(plotlist = vd$plotlist, ncol = 1)

app <- ExploreModelMatrix(sampleData = sampleData,
                          designFormula = ~ batch + condition)

if (interactive()) shiny::runApp(app)



## EJERCICIO

## Interpreta ResponseResistant.Treatmentpre del ejercicio 2. Puede ser útil
## tomar un screenshot (captura de pantalla) y anotarla con líneas de colores.
## Si haces eso, puedes incluir la imagen en tus notas.

## Interpretación matemática:
## E[y|treatment=pre,response=resistant,ind.n]-E[y|treatment=post,response=resistant,ind.n]

## Interpretación en palabras simples
## Se trata de la tasa promedio de cambio entre los niveles de expresión de en
## el mismo individuo antes (pre) y después del tratamiento (post) medido
## exclusivamente en personas con una respuesta resistente al tratamiento


## ¿Por qué es clave el 0 al inicio de la fórmula en el ejercicio 3?
## Porque con el 0 mantenemos nos permite no estimar el intercepto, e caso de no
## ponerlo se estima el intercepto, lo cual no es coveniente, por eso es
## importante poner el 0.

## Para proyectos incluir siempre violin plot combinado con boxplot y puntos
## para ubicar la distribución de los datos.

