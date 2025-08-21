#install.packages("igraph")
library(igraph)
library(tidyverse)
library(readxl)
library(dplyr)
library(Matrix)

datos <- read_excel("Derechos_concedidos3.xlsx")
names(datos)
df <- datos %>%
  rename(
    nombre_solicitante = "Nombre Solicitante",
  )

df_red <- df %>%
  select(nombre_solicitante, Cuenca) %>%
  filter(!is.na(nombre_solicitante), !is.na(Cuenca)) %>%
  distinct()

# --- Construcción bipartita eficiente desde matriz de incidencia dispersa ---
# Factores compactos para índices de filas (solicitantes) y columnas (cuencas)
f_solic <- factor(df_red$nombre_solicitante)
f_cuenc <- factor(df_red$Cuenca)

# Matriz de incidencia dispersa (1 si existe vínculo solicitante–cuenca)
M <- sparseMatrix(
  i = as.integer(f_solic),
  j = as.integer(f_cuenc),
  x = 1L,
  dims = c(nlevels(f_solic), nlevels(f_cuenc)),
  dimnames = list(levels(f_solic), levels(f_cuenc))
)

# Grafo bipartito; igraph setea V(g)$type: TRUE para filas (solicitantes), FALSE para columnas (cuencas)
g <- graph_from_incidence_matrix(M, weighted = TRUE)

# (Opcional) Evita plotear todo el grafo completo aquí; la visualización se hará en pasos siguientes.

#BIPARTITA
#Definir nodos
solicitantes <- unique(df_red$nombre_solicitante)
cuencas <- unique(df_red$Cuenca)
V(g)$type <- V(g)$name %in% solicitantes

# Verificar que la red es bipartita
is_bipartite(g)

# Visualización de la red bipartita
plot(g,
     vertex.color = ifelse(V(g)$type, "red", "lightgreen"),
     vertex.label = NA,
     vertex.size = ifelse(V(g)$type, 5, 8),
     edge.color = "gray70",
     layout = layout_as_bipartite(g),
     main = "Red Bipartita")


# --- Proyecciones eficientes usando álgebra de matrices dispersas ---
# Proyección solicitante–solicitante: W_s = M %*% t(M)
W_s <- tcrossprod(M)
diag(W_s) <- 0      # remover autoconexiones
W_s <- drop0(W_s)   # compactar ceros explícitos

# Proyección cuenca–cuenca: W_c = t(M) %*% M
W_c <- crossprod(M)
diag(W_c) <- 0
W_c <- drop0(W_c)

# Construir grafos ponderados a partir de las matrices de adyacencia dispersas
g_solicitantes <- graph_from_adjacency_matrix(W_s, mode = "undirected", weighted = TRUE, diag = FALSE)
g_cuencas      <- graph_from_adjacency_matrix(W_c, mode = "undirected", weighted = TRUE, diag = FALSE)

# (Opcional) Filtros ligeros para trabajar con subgrafos más manejables
# threshold <- 2
# g_solicitantes <- delete_edges(g_solicitantes, E(g_solicitantes)[weight < threshold])
# g_solicitantes <- delete_vertices(g_solicitantes, degree(g_solicitantes) == 0)
# g_cuencas      <- delete_edges(g_cuencas, E(g_cuencas)[weight < threshold])
# g_cuencas      <- delete_vertices(g_cuencas, degree(g_cuencas) == 0)

# (Opcional) Persistir objetos para reutilizar sin recalcular
# saveRDS(g_solicitantes, "g_solicitantes.rds")
# saveRDS(g_cuencas, "g_cuencas.rds")


