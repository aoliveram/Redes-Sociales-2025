rm(list=ls())
library(EGAnet)
library(qgraph)
library(psych)

set.seed(1234)

################################################################################
# Simulemos datos de síntomas
################################################################################

# Caso 1: Datos totalmente aleatorios
data_aleatoria <- data.frame(
  Anhedonia = rnorm(100, mean = 5, sd = 2),
  Insomnio = rnorm(100, mean = 4, sd = 2),
  Fatiga = rnorm(100, mean = 4, sd = 2),
  Irritabilidad = rnorm(100, mean = 6, sd = 2),
  Ansiedad = rnorm(100, mean = 7, sd = 2),
  Preocupación = rnorm(100, mean = 6, sd = 2),
  Nerviosismo = rnorm(100, mean = 7, sd = 2),
  Inquietud = rnorm(100, mean = 5, sd = 2)
)

# Caso 2: Datos intencionalmente correlacionados en dos dimensiones latentes: 
# Ansiedad y Depresión
n <- 200  # número de observaciones

# Generamos datos latentes
factor_ansiedad <- rnorm(n, mean = 0, sd = 1) # irritabilidad, ansiedad, preocupación, nerviosismo, inquietud
factor_depresion <- rnorm(n, mean = 0, sd = 1) # anhedonia, insomnio, fatiga

# Crear síntomas como combinaciones lineales de los factores latentes + ruido
data_correlacion <- data.frame(
  Anhedonia = 0.7 * factor_depresion + rnorm(n, mean = 0, sd = 0.3),
  Insomnio = 0.6 * factor_depresion + 0.4 * factor_ansiedad + rnorm(n, mean = 0, sd = 0.3),
  Fatiga = 0.8 * factor_depresion + rnorm(n, mean = 0, sd = 0.3),
  Irritabilidad = 0.7 * factor_ansiedad + rnorm(n, mean = 0, sd = 0.3),
  Ansiedad = 0.9 * factor_ansiedad + rnorm(n, mean = 0, sd = 0.3),
  Preocupacion = 0.8 * factor_ansiedad + rnorm(n, mean = 0, sd = 0.3),
  Nerviosismo = 0.85 * factor_ansiedad + rnorm(n, mean = 0, sd = 0.3),
  Inquietud = 0.75 * factor_ansiedad + rnorm(n, mean = 0, sd = 0.3)
)


################################################################################
# Normalizamos los datos
################################################################################

# Normalizamos para que todas las variables tengan media 0 y desviación estándar 1
data_aleatoria <- scale(data_aleatoria)
data_correlacion <- scale(data_correlacion)


################################################################################
# Aplicamos EGA
################################################################################

ega_aleatorio <- EGA(data_aleatoria)
ega_correlacion <- EGA(data_correlacion)

# Resultados
print(ega_aleatorio)
print(ega_correlacion)

################################################################################
# Observamos la red
################################################################################
# Visualizar la red utilizando qgraph
qgraph(ega_aleatorio$network, 
       layout = "spring", 
       labels = colnames(data), 
       label.cex = 0.8, 
       node.size = 7, 
       edge.color = "blue",
       title = "Red de síntomas: Ansiedad y Depresión")

# Visualizar la red utilizando qgraph
qgraph(ega_correlacion$network, 
       layout = "spring", 
       labels = colnames(data), 
       label.cex = 0.8, 
       node.size = 7, 
       edge.color = "blue",
       title = "Red de síntomas: Ansiedad y Depresión")


# Resumen de las comunidades
print(ega_aleatorio$dim.variables)
print(ega_correlacion$dim.variables)

################################################################################
# Usando otro paquete basado en técnicas de bootstrapping: bootnet
################################################################################

library(bootnet)
network <- estimateNetwork(data_correlacion, default = "pcor", threshold = "sig", alpha = 0.05)
plot(network, layout = "spring", gtheme = "colorblind")
centralityPlot(network, include = "all", scale = "raw0")

# bootstrapping
boot_resultado <- bootnet(network, nBoots = 1000, nCores = 4)

# Resumen y visualización de los resultados
summary(boot_resultado)
plot(boot_resultado)
