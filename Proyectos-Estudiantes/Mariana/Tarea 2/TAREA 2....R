library(readxl)
library(igraph)
library(dplyr)
library(bipartite)
library(ggplot2)
library(tidyr)

rm(list = ls())

setwd("C:/Users/DCCS3/Documents/REDES SOCIALES")

red = read_graph("red_vertebrados.graphml", format = "graphml")



######
### Comparacion con redes teoricas
######

# componente gigante
componentes <- clusters(red)
componentes

g <- which.max(componentes$csize) # identificamos la componente gigante
red2 <- induced.subgraph(red, which(componentes$membership == g)) # nos quedamos con el componente gigante

# Tamaño de la red
tam <- length(V(red2))
densidad <- edge_density(red2)

# Red erdos renyi
red_er <- erdos.renyi.game(tam, densidad)

# Red small world
red_sw <- watts.strogatz.game(dim = 1, size = tam, nei = 3, p = 0.1)

# Red barabasi
red_bar <- barabasi.game(n = tam, power = 1, m = 2, directed = F, algorithm = "psumtree")

#graficos
plot(red_bar, layout = layout_with_kk,
     vertex.label.font=7,
     vertex.label.cex = 1,
     vertex.label.color = "black",
     vertex.color = "lightblue",
     vertex.label.dist = 2)
nombres = 1:length(V(red2))


#### los 4 graficos para ver la forma grafica de cada red ###

plot(red2, layout = layout_with_kk,
     vertex.color = "lightblue",
     edge.color="darkgray",
     vertex.label = NA)

plot(red_bar, layout = layout_with_kk,
   vertex.color = "lightblue", edge.color="darkgray",   vertex.label = NA)

plot(red_sw, layout = layout_with_kk,
     vertex.color = "lightblue", edge.color="darkgray",   vertex.label = NA)

plot(red_er, layout = layout_with_kk,
     vertex.color = "lightblue",
     edge.color="darkgray",    vertex.label = NA)

# Configurar la ventana en 2 filas x 2 columnas
par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))  # mar ajusta márgenes

# Red real de vertebrados
plot(red2, layout = layout_with_kk,
     vertex.color = "lightblue",
     edge.color = "darkgray",
     vertex.label = NA,
     main = "Red vertebrados")

# Red Barabasi
plot(red_bar, layout = layout_with_kk,
     vertex.color = "lightblue",
     edge.color = "darkgray",
     vertex.label = NA,
     main = "Modelo Barabási")

# Red Small-World
plot(red_sw, layout = layout_with_kk,
     vertex.color = "lightblue",
     edge.color = "darkgray",
     vertex.label = NA,
     main = "Modelo Watts–Strogatz")

# Red Erdos-Renyi
plot(red_er, layout = layout_with_kk,
     vertex.color = "lightblue",
     edge.color = "darkgray",
     vertex.label = NA,
     main = "Modelo Erdős–Rényi")

### distancia media
distancia_red <- round(mean_distance(red2),3) 
distancia_red_er <- round(mean_distance(red_er),3)
distancia_red_sw <- round(mean_distance(red_sw),3)
distancia_red_bar <- round(mean_distance(red_bar),3)

### densidad
densidad_red <- edge_density(red2)
densidad_red_er <- edge_density(red_er)
densidad_sw <- edge_density(red_sw)
densidad_red_bar <- edge_density(red_bar)

### Tamaño (nodos y aristas) -----
n_vertices_red    <- igraph::vcount(red2)
n_edges_red       <- igraph::ecount(red2)

n_vertices_red_er <- igraph::vcount(red_er)
n_edges_red_er    <- igraph::ecount(red_er)

n_vertices_red_sw <- igraph::vcount(red_sw)
n_edges_red_sw    <- igraph::ecount(red_sw)

n_vertices_red_bar <- igraph::vcount(red_bar)
n_edges_red_bar    <- igraph::ecount(red_bar)

### Diámetro 
diametro_red     <- igraph::diameter(red2)
diametro_red_er  <- igraph::diameter(red_er)
diametro_red_sw  <- igraph::diameter(red_sw)
diametro_red_bar <- igraph::diameter(red_bar)

### centralidad closeness
centralidad_cercania_red <- igraph::closeness(red,normalized=T)
centralidad_cercania_red_er <- igraph::closeness(red_er,normalized=T)
centralidad_cercania_red_sw <- igraph::closeness(red_sw,normalized=T)
centralidad_cercania_red_bar <- igraph::closeness(red_bar,normalized=T)

### centralidad betweenness
centralidad_intermediacion_red <- igraph::betweenness(red, normalized = T)
centralidad_intermediacion_red_er <- igraph::betweenness(red_er, normalized = T)
centralidad_intermediacion_red_sw <- igraph::betweenness(red_sw, normalized = T)
centralidad_intermediacion_red_bar <- igraph::betweenness(red_bar, normalized = T)

### STARGAZER: métricas red original

# Armar data.frame para stargazer
tabla_red <- data.frame(
  N_vertices      = n_vertices_red,
  N_aristas       = n_edges_red,
  Diametro        = diametro_red,
  Distancia_media = distancia_red,
  Densidad        = densidad_red,
  Centralidad_closeness = centralidad_cercania_red,
  Centralidad_betweenness = centralidad_intermediacion_red
)


# Tabla descriptiva
stargazer(tabla_red, type = "text", title = "Resumen descriptivo", summary = TRUE)
# Mostrar en consola (texto)
stargazer::stargazer(
  tabla_red,
  type   = "text",
  summary = FALSE,
  title  = "Métricas básicas - Red original (red2)"
)


# tabla comparativa
tabla2 <- data.frame(
  Red             = c("Red vertebrados", "Erdos-Reny", "Watts Strogratz", "Barabasi"),
  N_vertices      = c(n_vertices_red, n_vertices_red_er, n_vertices_red_sw, n_vertices_red_bar),
  N_aristas       = c(n_edges_red, n_edges_red_er, n_edges_red_sw, n_edges_red_bar),
  Diametro        = c(diametro_red, diametro_red_er, diametro_red_sw, diametro_red_bar),
  Distancia_media = c(distancia_red, distancia_red_er, distancia_red_sw, distancia_red_bar),
  Densidad        = c(densidad_red, densidad_red_er, densidad_sw, densidad_red_bar),
  centr_betweenness = c(mean(centralidad_intermediacion_red), mean(centralidad_intermediacion_red_er),mean(centralidad_intermediacion_red_sw), mean(centralidad_intermediacion_red_bar)),
  centr_closeness = c(mean(centralidad_cercania_red),mean(centralidad_cercania_red_er), mean(centralidad_cercania_red_sw), mean(centralidad_cercania_red_bar))
)


tabla2
stargazer(tabla2, summary = F, type = "text", rownames = F )


######
## Implementación de simulaciones / bootstrapping
######

set.seed(123)  # Reproducibilidad

n_sim <- 1000

# Crear dataframes vacíos para almacenar resultados
sim_er <- data.frame(Distancia_media = numeric(n_sim),
                     Diametro = numeric(n_sim),
                     Densidad = numeric(n_sim),
                     Closeness = numeric(n_sim),
                     Betweenness = numeric(n_sim))

sim_sw <- data.frame(Distancia_media = numeric(n_sim),
                     Diametro = numeric(n_sim),
                     Densidad = numeric(n_sim),
                     Closeness = numeric(n_sim),
                     Betweenness = numeric(n_sim))

sim_ba <- data.frame(Distancia_media = numeric(n_sim),
                     Diametro = numeric(n_sim),
                     Densidad = numeric(n_sim),
                     Closeness = numeric(n_sim),
                     Betweenness = numeric(n_sim))

# Bucle de simulación
for (i in 1:n_sim) {
  # Red ER
  er <- igraph::erdos.renyi.game(n_vertices_red, p.or.m = densidad_red, type = "gnp", directed = FALSE)
  sim_er[i, ] <- c(
    Distancia_media = tryCatch(igraph::mean_distance(er), error = function(e) NA),
    Diametro = tryCatch(igraph::diameter(er), error = function(e) NA),
    Densidad = igraph::edge_density(er),
    Closeness = mean(igraph::closeness(er, normalized = TRUE), na.rm = TRUE),
    Betweenness = mean(igraph::betweenness(er, normalized = TRUE), na.rm = TRUE)
  )
  
  # Red SW
  sw <- igraph::watts.strogatz.game(1, size = n_vertices_red, nei = 3, p = 0.1)
  sim_sw[i, ] <- c(
    Distancia_media = tryCatch(igraph::mean_distance(sw), error = function(e) NA),
    Diametro = tryCatch(igraph::diameter(sw), error = function(e) NA),
    Densidad = igraph::edge_density(sw),
    Closeness = mean(igraph::closeness(sw, normalized = TRUE), na.rm = TRUE),
    Betweenness = mean(igraph::betweenness(sw, normalized = TRUE), na.rm = TRUE)
  )
  
  # Red BA
  ba <- igraph::barabasi.game(n = n_vertices_red, power = 1, m = 2, directed = FALSE, algorithm = "psumtree")
  sim_ba[i, ] <- c(
    Distancia_media = tryCatch(igraph::mean_distance(ba), error = function(e) NA),
    Diametro = tryCatch(igraph::diameter(ba), error = function(e) NA),
    Densidad = igraph::edge_density(ba),
    Closeness = mean(igraph::closeness(ba, normalized = TRUE), na.rm = TRUE),
    Betweenness = mean(igraph::betweenness(ba, normalized = TRUE), na.rm = TRUE)
  )
}

# Valor observado en la red real
valor_real <- data.frame(
  Distancia_media = distancia_red,
  Diametro = diametro_red,
  Densidad = densidad_red,
  Closeness = mean(centralidad_cercania_red),
  Betweenness = mean(centralidad_intermediacion_red)
)

# Función para obtener p-valor empírico (de una cola)
p_valor_empirico <- function(simulaciones, valor_obs) {
  # Eliminar NA
  simulaciones <- simulaciones[!is.na(simulaciones)]
  if(length(simulaciones) == 0) return(NA)  # Si todas las simulaciones son NA
  mean(simulaciones >= valor_obs, na.rm = TRUE)
}

# Calcular p-valores para cada métrica y cada modelo
metricas <- c("Distancia_media", "Diametro", "Densidad", "Closeness", "Betweenness")

comparacion <- data.frame(
  Metrica = metricas,
  Real = unlist(valor_real),
  ER_p = sapply(metricas, function(m) p_valor_empirico(sim_er[[m]], valor_real[[m]])),
  SW_p = sapply(metricas, function(m) p_valor_empirico(sim_sw[[m]], valor_real[[m]])),
  BA_p = sapply(metricas, function(m) p_valor_empirico(sim_ba[[m]], valor_real[[m]]))
)

print(comparacion)
stargazer(comparacion, summary = F, type = "text", rownames = F )


######
## Visualización de resultados

# Preparar data para ggplot
sim_er$modelo <- "Erdos-Renyi"
sim_sw$modelo <- "Small-World"
sim_ba$modelo <- "Barabasi"

simulaciones_total <- rbind(sim_er, sim_sw, sim_ba)

# Transformar a formato largo
sim_long <- pivot_longer(simulaciones_total, cols = -modelo, names_to = "metrica", values_to = "valor")

# Agregar valores reales para graficar
valores_reales_long <- valor_real |> 
  pivot_longer(cols = everything(), names_to = "metrica", values_to = "valor_real")

# Boxplot con líneas de valor real
ggplot(sim_long, aes(x = modelo, y = valor, fill = modelo)) +
  geom_boxplot(outlier.alpha = 0.1, alpha = 0.7) +
  facet_wrap(~metrica, scales = "free", ncol = 2) +
  geom_hline(data = valores_reales_long, aes(yintercept = valor_real), linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Comparación de métricas de red con modelos nulos (bootstraping)",
       y = "Valor", x = "Modelo nulo") +
  theme(legend.position = "none")



###################
## Red filtrada  ##
###################

setwd("C:/Users/DCCS3/Documents/REDES SOCIALES")

#Base de datos
datos = read_excel("datos_red_nido.xlsx", na="NA")
datos = as.data.frame(datos)

#Filtro para quitar UKN
datos = datos[ !(datos$ARBOL %in% "UKN"), ]

# Vectores con nombres de arboles y vertebrados
arboles <- sort(unique(datos$ARBOL))
aves <- sort(unique(datos$AVE))

# Crear matriz de adyacencia vacía
matriz <- matrix(0, nrow = length(aves), ncol = length(aves))
rownames(matriz) <- aves
colnames(matriz) <- aves

# Recorrer cada árbol
for (j in 1:length(arboles)) {
  arbol_actual <- arboles[j]
  
  # Subconjunto del árbol actual
  df_arbol <- subset(datos, ARBOL == arbol_actual)
  
  # Contar frecuencias por ave
  frec_aves <- table(df_arbol$AVE)
  
  # Filtrar solo aves con frecuencia >= 5
  aves_validas <- names(frec_aves[frec_aves >= 5])
  
  # Generar combinaciones si hay al menos 2 aves
  if (length(aves_validas) >= 2) {
    combinaciones <- combn(aves_validas, 2)
    
    for (k in 1:ncol(combinaciones)) {
      ave1 <- combinaciones[1, k]
      ave2 <- combinaciones[2, k]
      
      # asignar enlance
      matriz[ave1, ave2] <- 1
      matriz[ave2, ave1] <- 1
    }
  }
}

# Eliminar auto-conexiones
diag(matriz) <- 0

# Crear red no ponderada
red_filtrada <- igraph::graph_from_adjacency_matrix(matriz, mode = "undirected", diag = FALSE)

# Graficar red (sin labels)
plot(red_filtrada, 
     layout = layout_with_kk,
     vertex.label = NA,
     vertex.color = "lightblue",
     vertex.size = 10,
     edge.color = "gray")


# componente gigante
componentes2 <- clusters(red_filtrada)
g2 <- which.max(componentes2$csize) # identificamos la componente gigante
red_filtrada2 <- induced.subgraph(red, which(componentes2$membership == g2)) # nos quedamos con el componente gigante

plot(red_filtrada2, 
     layout = layout_with_kk,
     vertex.label = NA,
     vertex.color = "lightblue",
     vertex.size = 10,
     edge.color = "gray")
