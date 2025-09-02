library(igraph)
library(dplyr)
library(readr)
library(Matrix)
library(ggraph)
library(ggplot2)

#df <- read_csv("df_cite.csv")
df <- read_csv("Desktop/Doctorado/2025 - VI/Redes-Sociales-2025/Proyectos-Estudiantes/Jorge/Tarea 1/df_cite.csv")

autores_factor <- factor(df$Autor)
afiliaciones_factor <- factor(df$Afil)


M <- sparseMatrix(
  i = as.integer(autores_factor),
  j = as.integer(afiliaciones_factor),
  x = 1, 
  dims = c(nlevels(autores_factor), nlevels(afiliaciones_factor)),
  dimnames = list(levels(autores_factor), levels(afiliaciones_factor))
)


min_deg_autor <- 2L   
min_deg_afil <- 1L

rdeg <- Matrix::rowSums(M != 0) 
cdeg <- Matrix::colSums(M != 0) 

keep_rows  <- which(rdeg >= min_deg_autor)
keep_cols  <- which(cdeg >= min_deg_afil)

M_reducida <- M[keep_rows, keep_cols, drop = FALSE]
print(paste("Matriz original:", dim(M)[1], "autores x", dim(M)[2], "afiliaciones"))
print(paste("Matriz reducida:", dim(M_reducida)[1], "autores x", dim(M_reducida)[2], "afiliaciones"))

W_autores <- tcrossprod(M_reducida)
Matrix::diag(W_autores) <- 0L

W_autores@x[W_autores@x < umbral_peso] <- 0 
W_autores <- drop0(W_autores) 

W_afiliaciones <- crossprod(M_reducida) 
Matrix::diag(W_afiliaciones) <- 0L

W_afiliaciones@x[W_afiliaciones@x < umbral_peso] <- 0
W_afiliaciones <- drop0(W_afiliaciones)

g_autores <- graph_from_adjacency_matrix(W_autores, mode = "undirected", weighted = TRUE)
g_afiliaciones <- graph_from_adjacency_matrix(W_afiliaciones, mode = "undirected", weighted = TRUE)

gc_id_aut <- which.max(comp_aut$csize)
g_gc_aut  <- induced_subgraph(g_autores, which(comp_aut$membership == gc_id_aut))


V(g_gc_aut)$degree <- degree(g_gc_aut)
E(g_gc_aut)$width <- E(g_gc_aut)$weight

p_autores <- ggraph(g_gc_aut, layout = "fr") +
  geom_edge_link(aes(edge_width = width), alpha = 0.2, color = "gray60") +
  scale_edge_width_continuous(range = c(0.1, 1.5)) +
  geom_node_point(aes(size = degree), alpha = 0.7, color = "steelblue") +
  scale_size_continuous(range = c(0.5, 5)) +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Red de Coautoría")

print(p_autores)


comp_afil <- components(g_afiliaciones)
gc_id_afil <- which.max(comp_afil$csize)
g_gc_afil  <- induced_subgraph(g_afiliaciones, which(comp_afil$membership == gc_id_afil))

V(g_gc_afil)$degree <- degree(g_gc_afil)
E(g_gc_afil)$width <- E(g_gc_afil)$weight

p_afiliaciones <- ggraph(g_gc_afil, layout = "fr") +
  geom_edge_link(aes(edge_width = width), alpha = 0.4, color = "gray50") +
  scale_edge_width_continuous(range = c(0.2, 3)) +
  geom_node_point(aes(size = degree), color = "seagreen", alpha = 0.8) +
  geom_node_text(aes(label = ifelse(degree > quantile(degree, 0.8), name, '')), repel = TRUE, size = 3) +
  scale_size_continuous(range = c(1, 8)) +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Red de Colaboración Institucional")
  
print(p_afiliaciones)

comunidades_autores <- cluster_louvain(g_gc_aut)
V(g_gc_aut)$community <- as.factor(comunidades_autores$membership)
print(paste("Se encontraron", length(levels(V(g_gc_aut)$community)), "comunidades de autores."))

comunidades_afiliaciones <- cluster_louvain(g_gc_afil)
V(g_gc_afil)$community <- as.factor(comunidades_afiliaciones$membership)
print(paste("Se encontraron", length(levels(V(g_gc_afil)$community)), "comunidades de instituciones."))

p_autores_comunidades <- ggraph(g_gc_aut, layout = "fr") +
  geom_edge_link(alpha = 0.15) +
  geom_node_point(aes(size = degree, color = community), alpha = 0.8) +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Comunidades en la Red de Coautoría")
print(p_autores_comunidades)

p_afiliaciones_comunidades <- ggraph(g_gc_afil, layout = "fr") +
  geom_edge_link(aes(edge_width = width), alpha = 0.4, color = "gray50") +
  scale_edge_width_continuous(range = c(0.2, 3)) +
  geom_node_point(aes(size = degree, color = community), alpha = 0.8) +
  scale_size_continuous(range = c(1, 8)) +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Comunidades de Colaboración Institucional")
print(p_afiliaciones_comunidades)
