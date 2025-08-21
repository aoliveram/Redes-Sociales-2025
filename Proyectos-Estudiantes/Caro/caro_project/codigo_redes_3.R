library(igraph)
library(tidyr)
library(readxl)
library(dplyr)
library(Matrix)
library(ggplot2)
library(ggraph)

datos <- read_excel("Derechos_concedidos3.xlsx")  ## --> CAMBIAR ESTO !!
head(datos)

length(unique(datos$Región))

df <- datos %>%
  rename(
    nombre_solicitante = "Nombre Solicitante",
  )

df_red <- df %>%
  select(nombre_solicitante, Cuenca) %>%
  filter(!is.na(nombre_solicitante), !is.na(Cuenca)) %>%
  distinct()

# --- Construcción bipartita desde matriz de incidencia dispersa ---
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
g <- graph_from_biadjacency_matrix(M, weighted = TRUE)

# --- Gráfico: red bipartita completa (sin labels) ---
try({
  plot(
    g,
    vertex.color = ifelse(V(g)$type, "steelblue", "seagreen3"),
    vertex.label = NA,
    vertex.size = ifelse(V(g)$type, 2, 3),
    edge.color = "gray80",
    layout = layout_as_bipartite(g),
    main = "Red bipartita completa (solicitantes–cuencas)"
  )
}, silent = TRUE)

# Limpieza de nombres
rownames(M) <- NULL # nombres de cuencas se mantienen

# Filtrar por grado mínimo
min_deg_solic <- 2L   # solicitante debe aparecer en >= 2 cuencas
min_deg_cuenc <- 1L   # cuenca con al menos este número de solicitantes

rdeg <- Matrix::rowSums(M != 0)
cdeg <- Matrix::colSums(M != 0)

keep_rows  <- which(rdeg >= min_deg_solic)
keep_cols  <- which(cdeg >= min_deg_cuenc)

M_red <- M[keep_rows, keep_cols, drop = FALSE]

# --- Proyecciones sobre la matriz reducida ---

# Proyección solicitante–solicitante: W_s = M_red %*% t(M_red)
W_s <- tcrossprod(M_red)
Matrix::diag(W_s) <- 0L

w_threshold <- 2L  # solo se mantienen conexiones entre solicitantes que comparten al menos 2 cuencas
W_s@x[W_s@x < w_threshold] <- 0 # pone en cero los pesos menores al umbral
W_s <- drop0(W_s) # elimina entradas con peso cero

# Proyección cuenca–cuenca: W_c = t(M_red) %*% M_red
W_c <- crossprod(M_red)
Matrix::diag(W_c) <- 0L

W_c@x[W_c@x < w_threshold] <- 0
W_c <- drop0(W_c)

# --- Grafos ponderados a partir de proyecciones ---
g_solicitantes <- graph_from_adjacency_matrix(W_s, mode = "undirected", weighted = TRUE, diag = FALSE)
g_cuencas      <- graph_from_adjacency_matrix(W_c, mode = "undirected", weighted = TRUE, diag = FALSE)

# Para los solicitantes, asigna IDs cortos como nombres (evita strings largos)
V(g_solicitantes)$name <- as.character(seq_len(vcount(g_solicitantes)))

# --- Gráficos: proyecciones completas (sin labels) ---
# Proyección solicitantes–solicitantes

# Nos quedamos con el componente conexo más grande
comp_s <- components(g_solicitantes)
gc_id_s <- which.max(comp_s$csize)
g_gc_s  <- induced_subgraph(g_solicitantes, which(comp_s$membership == gc_id_s))

try({
  ggraph(g_gc_s, layout = "fr") +
    geom_edge_link(aes(width = weight), alpha = 0.15) +
    scale_edge_width(range = c(0.2, 1.5)) +
    geom_node_point(size = 1.2) +
    theme_void() +
    ggtitle("Proyección solicitantes–solicitantes (Componente Principal)")
  }, silent = TRUE)

# Proyección cuencas–cuencas

# Nos quedamos con el componente conexo más grande
comp_c <- components(g_cuencas)
gc_id_c <- which.max(comp_c$csize)
g_gc_c  <- induced_subgraph(g_cuencas, which(comp_c$membership == gc_id_c))

try({
  ggraph(g_gc_c, layout = "fr") +
    geom_edge_link(aes(width = weight), alpha = 0.25) +
    scale_edge_width(range = c(0.3, 2)) +
    geom_node_point(size = 2.2) +
    theme_void() +
    ggtitle("Proyección cuencas–cuencas (Componente Principal)")
  }, silent = TRUE)
