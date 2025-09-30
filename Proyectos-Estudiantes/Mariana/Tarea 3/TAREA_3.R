library(readxl)
library(igraph)
library(bipartite)
library(stargazer)

setwd("C:/Users/DCCS3/Documents/REDES SOCIALES/tarea 3/")

####
### para fines prácticos y en vista de que conforman la mayoria de las especies,
### llamaré "AVES" a todo el conjunto de vertabrados en este codigo  ###
####

rm(list = ls())

#Base de datos
datos = read_excel("DATA_DECAY.xlsx", na="NA")
datos = datos[ , c(3,7) ]
datos = na.omit(datos)
sanos = datos[datos$DECAY %in% c("2A", "2B", "2C"),  ]
muertos = datos[datos$DECAY %in% c("3","4","5","6","7","8"),  ]

## --- DESCRIPTIVO: todos los DECAY + Sanos(vivos) vs Muertos

# Orden explícito de categorías DECAY
dec_levels <- c("2A","2B","2C","3","4","5","6","7","8")

# Conteo por categoría (incluye 0 si falta alguna)
tab_all <- table(factor(datos$DECAY, levels = dec_levels))

# Agregación Sanos (2A–2c) vs Muertos (3–8)
sanos_levels   <- c("2A","2B","2C")
muertos_levels <- c("3","4","5","6","7","8")
counts_two <- c(
  Sanos   = sum(tab_all[sanos_levels],   na.rm = TRUE),
  Muertos = sum(tab_all[muertos_levels], na.rm = TRUE)
)

# Graficar lado a lado
op <- par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

# 1) Barras por cada DECAY
barplot(
  tab_all,
  main = "Frecuencia por DECAY",
  xlab = "Categoría DECAY",
  ylab = "N° de registros",
  las = 1,
  cex.names = 1.5,   # tamaño etiquetas X (categorías DECAY)
  cex.axis  = 1.2,   # tamaño números eje Y
  cex.lab   = 1.4,   # tamaño textos ejes
  cex.main  = 1.5    # tamaño título
)

# 2) Barras Sanos vs Muertos
barplot(
  counts_two,
  main = "Árboles sanos vs Árboles muertos",
  ylab = "N° de registros",
  las = 1,
  names.arg = c("Sanos (2A– 2C)", "Muertos (3–8)"),
  col = c("cyan4", "bisque3"),
  cex.names = 1.5,   # etiquetas bajo cada barra
  cex.axis  = 1.2,   # números del eje Y
  cex.lab   = 1.4,   # ejes
  cex.main  = 1.5    # título
)
par(op)

sanos
muertos
## --- Matriz ponderada AVE × DECAY ---
# (peso = frecuencia de interacciones AVE–DECAY)
dec_levels <- c("2A","2B","2C","3","4","5","6","7","8")   # orden deseado

aves_levels <- unique(datos$AVE)
aves_levels
tab <- with(datos, table(factor(AVE,  levels = aves_levels),
                         factor(DECAY, levels = dec_levels)))
M <- as.matrix(tab)

## --- Orden anidado ---
### orden decreciente anidado, las más conectadas estan en la parte superior.
###
M_nested <- sortweb(M, sort.order = "dec")   

## --- Plot bipartito con pesos ---
op <- par(mar = c(6,6,7 ,5 ))  # márgenes más amplios para etiquetas
plotweb(
  M,
  method   = "normal",
  text.rot = 90,
  col.low  = "cyan",
  col.high = "coral4",
  ybig     = 1.8,
  low.lablength  = 6,
  high.lablength = 6,
  labsize = 1.2   # controla tamaño global de etiquetas
)
par(op)


#RED SANOS

# Vectores únicos
decay <- sort(unique(sanos$DECAY))
aves <- sort(unique(sanos$AVE))

# Crear matriz de adyacencia vacía
matriz <- matrix(0, nrow = length(aves), ncol = length(aves))
rownames(matriz) <- aves
colnames(matriz) <- aves

# Recorrer cada árbol
for (j in 1:length(decay)) {
  decay_actual <- decay[j]
  
  # Subconjunto del árbol actual
  df_decay <- subset(sanos, DECAY == decay_actual)
  
  # Contar frecuencias por ave
  frec_aves <- table(df_decay$AVE)
  
  # Filtrar solo aves con frecuencia >= 10
  aves_validas <- names(frec_aves[frec_aves >= 3])
  
  # Generar combinaciones si hay al menos 2 aves
  if (length(aves_validas) >= 2) {
    combinaciones <- combn(aves_validas, 2)
    
    for (k in 1:ncol(combinaciones)) {
      ave1 <- combinaciones[1, k]
      ave2 <- combinaciones[2, k]
      
      # Marcar relación binaria (sin peso)
      matriz[ave1, ave2] <- 1
      matriz[ave2, ave1] <- 1
    }
  }
}

# Eliminar auto-conexiones
diag(matriz) <- 0

# Crear red no ponderada
red_sanos <- igraph::graph_from_adjacency_matrix(matriz, mode = "undirected", diag = FALSE)

#RED MUERTOS

# Vectores únicos
decay <- sort(unique(muertos$DECAY))
aves <- sort(unique(muertos$AVE))

# Crear matriz de adyacencia vacía
matriz <- matrix(0, nrow = length(aves), ncol = length(aves))
rownames(matriz) <- aves
colnames(matriz) <- aves

# Recorrer cada árbol
for (j in 1:length(decay)) {
  decay_actual <- decay[j]
  
  # Subconjunto del árbol actual
  df_decay <- subset(muertos, DECAY == decay_actual)
  
  # Contar frecuencias por ave
  frec_aves <- table(df_decay$AVE)
  
  # Filtrar solo aves con frecuencia >= 10
  aves_validas <- names(frec_aves[frec_aves >= 3])
  
  # Generar combinaciones si hay al menos 2 aves
  if (length(aves_validas) >= 2) {
    combinaciones <- combn(aves_validas, 2)
    
    for (k in 1:ncol(combinaciones)) {
      ave1 <- combinaciones[1, k]
      ave2 <- combinaciones[2, k]
      
      # Marcar relación binaria (sin peso)
      matriz[ave1, ave2] <- 1
      matriz[ave2, ave1] <- 1
    }
  }
}

# Eliminar auto-conexiones
diag(matriz) <- 0

# Crear red no ponderada
red_muertos <- igraph::graph_from_adjacency_matrix(matriz, mode = "undirected", diag = FALSE)


# Configurar la ventana en 2 filas x 2 columnas
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))  # mar ajusta márgenes

##para los sanos

excavadores <- c("PYGALB", "COLPIT", "VENLIG", "CAMMAG")
colores <- rep("cyan3", igraph::vcount(red_sanos))
indices <- which(igraph::V(red_sanos)$name %in% excavadores)
colores[indices] <- "coral"

##para los muertos en pie

excavadores <- c("PYGALB", "COLPIT", "VENLIG", "CAMMAG")
colores2 <- rep("cyan3", igraph::vcount(red_muertos))
indices2 <- which(igraph::V(red_muertos)$name %in% excavadores)
colores2 [indices2] <- "coral"

# Graficar red sanos
plot(red_sanos, 
     layout = layout_with_kk,
     vertex.color = colores,
     vertex.size = 10,
     edge.color = "gray",
     main = " Árboles sanos")

# Graficar red muertos
plot(red_muertos, 
     layout = layout_with_kk,
     vertex.color = colores2,
     vertex.size = 10,
     edge.color = "gray",
  main = "Árboles muertos")


# Redes sin los nodos sueltos

nodos_aislados <- V(red_sanos)[igraph::degree(red_sanos) == 0]
red_sanos2 <- delete_vertices(red_sanos, nodos_aislados)

nodos_aislados2 <- V(red_muertos)[igraph::degree(red_muertos) == 0]
red_muertos2 <- delete_vertices(red_muertos, nodos_aislados2)

# Graficar red sanos
plot(red_sanos2, 
     layout = layout_with_kk,
     vertex.size = 10,
     vertex.color = "cyan3",
     edge.color = "gray",
     main = " Árboles sanos")

# Graficar red muertos
plot(red_muertos2, 
     layout = layout_with_kk,
     vertex.size = 10,
     vertex.color = "cyan3",
     edge.color = "gray",
     main = "Árboles muertos")

# Crear una función para extraer métricas clave
obtener_metricas <- function(g) {
  comps <- igraph::components(g)
  data.frame(
    nodos = igraph::vcount(g),
    enlaces = igraph::ecount(g),
    densidad = igraph::edge_density(g),
    grado_medio = mean(igraph::degree(g)),
    transitividad = igraph::transitivity(g, type = "global"),
    modularidad = igraph::modularity(igraph::cluster_fast_greedy(g)),
    diametro = ifelse(comps$no > 1, NA, igraph::diameter(g)),
    componentes = comps$no
  )
}

# Calcular métricas para ambas redes
m_sanos2 <- obtener_metricas(red_sanos2)
m_muertos2 <- obtener_metricas(red_muertos2)

# Muestra de resultados
# Comparar lado a lado
tabla_resultados = rbind(Sanos = m_sanos2, Muertos = m_muertos2)

stargazer::stargazer(
  tabla_resultados,
  type   = "text",
  summary = FALSE,
  title  = "Métricas básicas"
)


############
####### Diferencias en el grado medio
##########

comparar_metrica <- function(g1, g2, metrica_fun, n = 1000) {
  m1 <- metrica_fun(g1)
  m2 <- metrica_fun(g2)
  obs_diff <- abs(m1 - m2)
  
  g_total <- igraph::disjoint_union(g1, g2)
  nodos_totales <- igraph::vcount(g_total)
  
  difs_perm <- numeric(n)
  for (i in 1:n) {
    nodos_perm <- sample(igraph::V(g_total))
    g1_perm <- igraph::induced_subgraph(g_total, nodos_perm[1:igraph::vcount(g1)])
    g2_perm <- igraph::induced_subgraph(g_total, nodos_perm[(igraph::vcount(g1) + 1):nodos_totales])
    
    d1 <- tryCatch(metrica_fun(g1_perm), error = function(e) NA)
    d2 <- tryCatch(metrica_fun(g2_perm), error = function(e) NA)
    
    difs_perm[i] <- abs(d1 - d2)
  }
  
  p_val <- mean(difs_perm >= obs_diff, na.rm = TRUE)
  
  return(list(
    obs_diff = obs_diff,
    p_val = p_val,
    distrib = difs_perm
  ))
}

# Aplicar al grado medio
set.seed(123)
resultado_grado <- comparar_metrica(
  red_sanos2,
  red_muertos2,
  metrica_fun = function(g) mean(igraph::degree(g)),
  n = 1000
)

print(resultado_grado)


#######
### comparación con redes teoricas
#####

comparar_modelos <- function(red_gigante, nombre_red) {
  tam <- length(igraph::V(red_gigante))
  dens <- igraph::edge_density(red_gigante)
  
  # Generar modelos teóricos
  red_er <- igraph::sample_gnp(n = tam, p = dens)
  red_sw <- igraph::sample_smallworld(dim = 1, size = tam, nei = 2, p = 0.1)
  red_bar <- igraph::sample_pa(n = tam, m = round(dens * tam), directed = FALSE)
  
  # Función auxiliar para extraer métricas
  metricas <- function(red) {
    c(
      Densidad = round(igraph::edge_density(red), 3),
      Diametro = igraph::diameter(red),
      Clustering = round(mean(igraph::transitivity(red, type = "local", isolates = "zero")), 3)
    )
  }
  
  ###
  #### Crear tabla con resultados
  tabla_resultados <- rbind(
    Real            = metricas(red_gigante),
    Erdos_Renyi     = metricas(red_er),
    Small_World     = metricas(red_sw),
    Barabasi_Albert = metricas(red_bar)
  )
  
  # Convertir a data.frame y agregar columna para identificar modelo
  df_resultados <- as.data.frame(tabla_resultados)
  df_resultados$modelo <- rownames(df_resultados)
  rownames(df_resultados) <- NULL
  # Reordenar columnas para que 'modelo' quede primero
  df_resultados <- df_resultados[, c("modelo", names(df_resultados)[-ncol(df_resultados)])]
  
  return(df_resultados)
}

m_sanos2 <- comparar_modelos(red_sanos_gigante, "red_sanos2")
m_muertos2 <- comparar_modelos(red_muertos_gigante, "red_muertos2")

# Combinar resultados de ambas redes para mostrar
tabla_combinada <- rbind(
  Sanos2 = m_sanos2,
  Muertos2 = m_muertos2
)

# Para evitar duplicar la columna 'modelo' (revisar estructura)
tabla_combinada$modelo <- paste0(rownames(tabla_combinada), "_", tabla_combinada$modelo)
rownames(tabla_combinada) <- NULL

# Mostrar con stargazer
stargazer::stargazer(
  tabla_combinada,
  type    = "text",
  summary = FALSE,
  title   = "Métricas básicas para redes Sanos2 y Muertos2",
  rownames = FALSE
)