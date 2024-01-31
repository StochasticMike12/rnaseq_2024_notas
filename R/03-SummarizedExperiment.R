################################################################################
#                                                                              #
#     3. Objetos de Bioconductor para datos de expresión        31/01/2024     #
#                                                                              #
################################################################################

## Lets build our first SummarizedExperiment object
library("SummarizedExperiment")
## ?SummarizedExperiment

## De los ejemplos en la ayuda oficial

## Creamos los datos para nuestro objeto de tipo SummarizedExperiment
## para 200 genes a lo largo de 6 muestras
nrows <- 200
ncols <- 6
## Números al azar de cuentas
set.seed(20210223)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
## Información de nuestros genes
rowRanges <- GRanges(
  rep(c("chr1", "chr2"), c(50, 150)),
  IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
  strand = sample(c("+", "-"), 200, TRUE),
  feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))
## Información de nuestras muestras
colData <- DataFrame(
  Treatment = rep(c("ChIP", "Input"), 3),
  row.names = LETTERS[1:6]
)
## Juntamos ahora toda la información en un solo objeto de R
rse <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  rowRanges = rowRanges,
  colData = colData
)

## Exploremos el objeto resultante
rse



## 3.2. EJERCICIO

## Comando 1
rse[1:2, ]
## El comando nos muestra el resumen de los primeros 6 genes en todas las
## muestras.


## Comando 2
rse[, c("A", "D", "F")]
## El comando nos muestra el resumen de todos los genes en las meuestras A, D y
## F.


## Función para eliminar objetos
rm()

## Función para mostrar los objetos existentes
ls()



## 3.4 EJERCICIO CON SPATIALLIBD

## Descarguemos unos datos de spatialLIBD
sce_layer <- spatialLIBD::fetch_data("sce_layer")

sce_layer

## Revisemos el tamaño de este objeto
lobstr::obj_size(sce_layer)

iSEE::iSEE(sce_layer)


## Explora en con un heatmap la expresión de los genes MOBP, MBP y PCP4. Si
## hacemos un clustering (agrupamos los genes), ¿cúales genes se parecen más?
## Los patrones de expresión de los genes MOBP y MBP son más parecidos.

## ¿En qué capas se expresan más los genes MOBP y MBP?
## Se expresan más en la capa WM



