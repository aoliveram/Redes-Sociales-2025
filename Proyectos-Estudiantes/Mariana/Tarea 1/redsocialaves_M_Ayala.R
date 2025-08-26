library(readxl)
library(igraph)
library(bipartite)

rm(list = ls())

setwd("C:/Users/DCCS3/Documents/REDES SOCIALES")

#Base de datos
datos = read_excel("datos_red_nido.xlsx", na="NA")
datos = as.data.frame(datos)

# Filtrar para que sean solo aves
# filtronoaves = c("ANAFLA", "MURCI","DROGLI","LIOCHI","LIOTEN", "LEPEAG","LASVAR","LEOGUI", "LIOPIC", "RATRAT")
# datos = datos[ !(datos$AVE %in% filtronoaves), ]


##########################
### RED DE VERTEBRADOS ###
##########################

#Vector arboles contiene todos los posibles nombres de especies de arboles de la base de datos
arboles = unique(datos$ARBOL)

#Vector aves contiene todos los posibles nombres de especies de aves de la base de datos
aves = unique(datos$AVE)

# Matriz vacia para luego indicar conexiones
matriz = matrix(data = 0, nrow = length(aves), ncol = length(aves))
colnames(matriz) = aves
rownames(matriz) = aves

for (j in 1:length(arboles)) {
  # Subconjunto de aves asociadas a este árbol
  df = datos[ datos$ARBOL == arboles[j], ]
  aves_en_arbol = unique(df$AVE)
  
  # Todas las combinaciones de aves de a pares
  if (length(aves_en_arbol) >= 2) {
    combinaciones = combn(aves_en_arbol, 2)
    
    for (k in 1:ncol(combinaciones)) {
      ave1 = combinaciones[1, k]
      ave2 = combinaciones[2, k]
      
      # Marcar relación en la matriz (simétrica)
      matriz[ave1, ave2] = 1
      matriz[ave2, ave1] = 1
    }
  }
}

# Diagonal a 0 para no contar relaciones consigo mismas
diag(matriz) = 0

red <- graph_from_adjacency_matrix(matriz, mode = "undirected", diag = FALSE)

plot(red, layout = layout_with_kk,
     vertex.label.font=7,
     vertex.label.cex = 1,
     vertex.label.color = "black",
     vertex.color = "lightblue",
     vertex.size = 10,
     vertex.label.dist = 2)

################################################
### Red con nodos de tamaño proporcional #######
### al número de datos del mismo vertebrado ####
################################################

# Calcular el número de veces que cada especie de ave aparece en los datos
frecuencias_aves = table(datos$AVE)

# Normalizar las frecuencias para obtener un valor entre 1 y 50
tamanos_aves = (frecuencias_aves / sum(frecuencias_aves)) * 30  # 50 es el tamaño máximo de los vértices

# Obtener el tamaño de los nodos según la frecuencia normalizada
vertex_sizes = tamanos_aves[rownames(matriz)]

# Graficar la red
plot(red, layout = layout_with_kk,
     vertex.label.font=7,
     vertex.label.cex = 1,
     vertex.label.color = "black",
     vertex.color = "lightblue",
     vertex.size = vertex_sizes,  # Usar los tamaños calculados
     vertex.label.dist = 2)


########################
#### RED BIPARTITA  ####
########################

# Crear una matriz de adyacencia con filas para vertebrados y columnas para árboles
matriz_bipartita = matrix(0, nrow = length(aves), ncol = length(arboles))
rownames(matriz_bipartita) = aves
colnames(matriz_bipartita) = arboles

# Llenar la matriz con 1 si hay relación entre ave y árbol
for (j in 1:length(arboles)) {
  # Subconjunto de aves asociadas a este árbol
  df = datos[datos$ARBOL == arboles[j], ]
  aves_en_arbol = unique(df$AVE)
  
  for (ave in aves_en_arbol) {
    # Marcar la relación en la matriz (1 si ave está asociada al árbol)
    matriz_bipartita[ave, arboles[j]] = 1
  }
}

# Grafo bipartito: se indica que las dos partes son aves y árboles
red_bipartita <- graph_from_incidence_matrix(matriz_bipartita)

# Graficar la red bipartita
plot(red_bipartita,
     layout = layout_as_bipartite,
     vertex.label.font = 7,
     vertex.label.cex = 0.8,
     vertex.label.color = "black",
     vertex.color = c(rep("lightblue", length(aves)), rep("lightgreen", length(arboles))),
     vertex.size = 5,
     vertex.label.dist = 2,
     vertex.label.degree = -pi/2,
     asp = 0)

###########################################
## RED BIPARTITA VERTEBRADO-ARBOL IGRAPH ##
###########################################

ab_ave <- table(datos$AVE)  # Cuenta cuántas veces aparece cada ave
ab_arbol <- table(datos$ARBOL)  # Cuenta cuántas veces aparece cada árbol

# Obtener el vector de conteos para aves y árboles
ab_vec <- ifelse(V(red_bipartita)$type, ab_arbol[V(red_bipartita)$name], ab_ave[V(red_bipartita)$name])
ab_vec[is.na(ab_vec)] <- 0

# Función de reescalado manual
rescale <- function(x, to = c(8, 28)) {
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) return(rep(mean(to), length(x)))
  (x - r[1]) / diff(r) * diff(to) + to[1]
}

# Asignar tamaño de los vértices según los conteos
V(red_bipartita)$size  <- ifelse(V(red_bipartita)$type, rescale(ab_vec, c(14, 34)), rescale(ab_vec, c(8, 26)))
V(red_bipartita)$shape <- ifelse(V(red_bipartita)$type, "square", "circle")
V(red_bipartita)$color <- ifelse(V(red_bipartita)$type, "lightgreen", "lightblue")
V(red_bipartita)$frame.color <- NA
V(red_bipartita)$label <- V(red_bipartita)$name
V(red_bipartita)$label.cex <- rescale(V(red_bipartita)$size, to = c(0.7, 1.5))
V(red_bipartita)$label.color <- "black"
V(red_bipartita)$label.font  <- 2

# Layout con Fruchterman-Reingold inicial, luego ajustando posiciones
set.seed(42)
coords0 <- layout_as_bipartite(red_bipartita, types = V(red_bipartita)$type)  # aves/árboles separados
coords_kk <- layout_with_kk(red_bipartita, coords = coords0, weights = E(red_bipartita)$weight)

# Separar más aves abajo / árboles arriba
coords_kk[, 2] <- ifelse(V(red_bipartita)$type, coords_kk[, 2] + 1.8, coords_kk[, 2] - 1.8)
coords_kk <- coords_kk * 1.3
rownames(coords_kk) <- V(red_bipartita)$name

# Graficar
op <- par(mar = c(5, 3, 4, 3))  # Ajustar márgenes de la gráfica
plot(red_bipartita, layout = coords_kk)


###########################
## BIPARTITA CON PLOTWEB ##
###########################

plotweb(sortweb(matriz_bipartita, sort.order="dec"), text.rot = 90, method="normal", col.high = "lightgreen", col.low = "lightblue") 
