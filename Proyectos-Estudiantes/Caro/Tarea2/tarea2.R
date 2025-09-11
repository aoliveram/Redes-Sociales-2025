library(igraph)
library(tidyr)
library(readxl)
library(dplyr)
library(Matrix)
library(ggplot2)
library(ggraph)

setwd("D:/OneDrive/DCCS/1er año/Redes Sociales/Tarea1")
datos <- read_excel("Derechos_concedidos4.xlsx")  
names(datos)
df <- datos %>%
  rename(
    nombre_solicitante = "Nombre Solicitante",
    caudal_anual = "Caudal \r\nAnual\r\nProm \r\n"
  )
df <- df %>%
  arrange(desc("caudal_anual")) %>%  
  slice_head(n = 10000)     
df_red <- df %>%
  select(nombre_solicitante, Cuenca, caudal_anual, Región) %>%
  filter(!is.na(nombre_solicitante), !is.na(Cuenca)) %>%
  distinct()

f_solic <- factor(df_red$nombre_solicitante)
f_cuenc <- factor(df_red$Cuenca)

M <- sparseMatrix(
  i = as.integer(f_solic),
  j = as.integer(f_cuenc),
  x = 1L,
  dims = c(nlevels(f_solic), nlevels(f_cuenc)),
  dimnames = list(levels(f_solic), levels(f_cuenc))
)

g <- graph_from_incidence_matrix(M, weighted = TRUE)

# --- Gráfico: red bipartita
layout_bip <- layout_as_bipartite(g, hgap = 3, vgap = 1)
col_solicitantes <- adjustcolor("steelblue", alpha = 0.7)
col_cuencas <- adjustcolor("darkorange2", alpha = 0.8)
grados <- degree(g)
tam_vertices <- ifelse(V(g)$type, 
                       pmax(0.5, pmin(4, sqrt(grados/10))), 
                       pmax(1, pmin(6, sqrt(grados/5))))

if("caudal_anual" %in% names(df_red)) {
  df_edges <- df_red %>%
    mutate(
      from = as.integer(f_solic),
      to = as.integer(f_cuenc) + nlevels(f_solic)
    ) %>%
    select(from, to, caudal_anual)
  edge_list <- get.edgelist(g)
  for(i in 1:nrow(edge_list)) {
    E(g)[i]$weight <- 1  
  }
}

plot(g,
     vertex.color = ifelse(V(g)$type, col_solicitantes, col_cuencas),
     vertex.label = NA,
     vertex.size = tam_vertices,
     vertex.frame.color = adjustcolor("white", alpha = 0.8),
     vertex.frame.width = 0.5,
     edge.color = adjustcolor("gray40", alpha = 0.25),
     edge.width = 0.3,
     layout = layout_bip,
     cex.main = 1.1)

rownames(M) <- NULL
M_red <- M

# --- Proyecciones sobre la matriz reducida ---
W_s <- tcrossprod(M_red)
Matrix::diag(W_s) <- 0L

w_threshold <- 2L  
W_s@x[W_s@x < w_threshold] <- 0 
W_s <- drop0(W_s) 

W_c <- crossprod(M_red)
Matrix::diag(W_c) <- 0L

W_c@x[W_c@x < w_threshold] <- 0
W_c <- drop0(W_c)

# --- Grafos ponderados a partir de proyecciones ---
g_solicitantes <- graph_from_adjacency_matrix(W_s, mode = "undirected", weighted = TRUE, diag = FALSE)
g_cuencas      <- graph_from_adjacency_matrix(W_c, mode = "undirected", weighted = TRUE, diag = FALSE)

V(g_solicitantes)$name <- as.character(seq_len(vcount(g_solicitantes)))

# --- Gráficos: proyecciones completas ---
# Proyección solicitantes–solicitantes
# Nos quedamos con el componente conexo más grande
comp_s <- components(g_solicitantes)
gc_id_s <- which.max(comp_s$csize)
g_gc_s  <- induced_subgraph(g_solicitantes, which(comp_s$membership == gc_id_s))
grados_nodos <- degree(g_gc_s)
min_size <- 0.8   # tamaño mínimo de nodo
max_size <- 4.0   # tamaño máximo de nodo
grados_normalizados <- scales::rescale(grados_nodos, to = c(min_size, max_size))
pesos_aristas <- E(g_gc_s)$weight
min_width <- 0.05  
max_width <- 2.0  
pesos_normalizados <- scales::rescale(pesos_aristas, to = c(min_width, max_width))

V(g_gc_s)$size <- grados_normalizados
V(g_gc_s)$degree <- grados_nodos  
E(g_gc_s)$width <- pesos_normalizados

try({
  p1 <- ggraph(g_gc_s, layout = "fr") +
    geom_edge_link(aes(edge_width = weight), alpha = 0.4, color = "gray60") +
    scale_edge_width_continuous(
      range = c(min_width, max_width)) +
    geom_node_point(aes(size = degree), alpha = 0.8, color = "steelblue") +
    scale_size_continuous(
      range = c(min_size, max_size)) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none") +
    ggtitle("Proyección Solicitantes–Solicitantes")
  print(p1)
}, silent = TRUE)


deg_gc <- degree(g_gc_s)
umbral_alto <- 1500  
nodos_filtrados <- V(g_gc_s)[deg_gc > 0 & deg_gc <= umbral_alto]
g_gc_s_filt <- induced_subgraph(g_gc_s, vids = nodos_filtrados)

###METRICAS
g <- simplify(g_gc_s_filt)  

#Tamaño red
num_nodos <- vcount(g)
num_aristas <- ecount(g)
#Diámetro
diametro <- diameter(g, directed = FALSE, unconnected = FALSE)
#Centralidades
#Centralidad de grado
cent_grado <- degree(g, mode = "all", normalized = TRUE)
#Centralidad de intermediación
cent_bet <- betweenness(g, directed = FALSE, normalized = TRUE)
#Centralidad de cercanía
cent_close <- closeness(g, normalized = TRUE)
#Centralidad eigenvector
cent_eigen <- eigen_centrality(g, directed = FALSE)$vector
#Distribución de grado ---
deg_dist <- degree_distribution(g, cumulative = FALSE)
#Densidad ---
densidad <- edge_density(g)
#Clusterización ---
#clustering local
clust_local <- transitivity(g, type = "local", isolates = "zero")
clust_media <- mean(clust_local, na.rm = TRUE)
#Clusterización global
clust_global <- transitivity(g, type = "global")
#Longitud de camino medio
long_camino_medio <- mean_distance(g, directed = FALSE, unconnected = FALSE)
long_camino_medio

#TABLA RESUMEN
tabla_resumen <- tibble(
  Metrica = c("Número de nodos", 
              "Número de aristas", 
              "Diámetro", 
              "Densidad", 
              "Clusterización global"),
  Valor = c(num_nodos, 
            num_aristas, 
            diametro, 
            round(densidad, 4), 
            round(clust_global, 4))
)
print(tabla_resumen)

#TABLA POR NODOS (TOP 10)
tabla_nodos <- tibble(
  Nodo = V(g)$name,
  Grado = cent_grado,
  Betweenness = cent_bet,
  Closeness = cent_close,
  Eigenvector = cent_eigen,
  Cluster_local = clust_local
) %>%
  arrange(desc(Grado)) %>%
  slice_head(n = 10)

print(tabla_nodos)

##CENTRALIDADES (visualizacion)
V(g)$cent_grado  <- cent_grado
V(g)$cent_bet    <- cent_bet
V(g)$cent_close  <- cent_close
V(g)$cent_eigen  <- cent_eigen

layout_fr <- create_layout(g, layout = "fr")
#centralidad de grado
try({
  p3 <- ggraph(layout_fr) +
    geom_edge_link(alpha = 0.4, color = "gray60") +
    geom_node_point(aes(size = degree, color = cent_grado), alpha = 0.8) +
    scale_size_continuous(
      range = c(min_size, max_size)) +
    scale_color_viridis_c(option = "plasma") + 
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right") +
    ggtitle("Centralidad de grado")
  print(p3)
}, silent = TRUE)
#centralidad betweness
try({
  p4 <- ggraph(layout_fr) +
    geom_edge_link(alpha = 0.4, color = "gray60") +
    geom_node_point(aes(size = degree, color = cent_bet), alpha = 0.8) +
    scale_size_continuous(
      range = c(min_size, max_size)) +
    scale_color_viridis_c(option = "plasma") + 
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right") +
    ggtitle("Centralidad Betweness")
  print(p4)
}, silent = TRUE)

#centralidad de eigenvector
try({
  p5 <- ggraph(layout_fr) +
    geom_edge_link(alpha = 0.4, color = "gray60") +
    geom_node_point(aes(size = degree, color = cent_eigen), alpha = 0.8) +
    scale_size_continuous(
      range = c(min_size, max_size)) +
    scale_color_viridis_c(option = "plasma") + 
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right") +
    ggtitle("Centralidad Eigenvector")
  print(p5)
}, silent = TRUE)

#Distribucion de grado (graficos)
# Vector de grados
deg <- degree(g)
# Gráfico de densidad
ggplot(data.frame(deg), aes(x = deg)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribución de Grado",
       x = "Grado (k)",
       y = "Frecuencia")

ggplot(data.frame(deg), aes(x = deg)) +
  geom_density(color = "black", fill = "steelblue", alpha = 0.4, adjust = 1.5, size = 0.5) +
  theme_minimal() +
  labs(
    title = "Densidad de Distribución de Grado",
    x = "Grado (k)",
    y = "Densidad (P(k))"
  )

#log-log
ggplot(data.frame(deg), aes(x = deg)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  theme_minimal() +
  scale_x_log10() + 
  scale_y_log10() +
  labs(title = "Distribución de Grado (log-log)",
       x = "Grado_k (log)",
       y = "Frecuencia (log)")

d <- density(deg)
df_d <- data.frame(x = d$x, y = d$y)
ggplot(df_d, aes(x = x, y = y)) +
  geom_line(color = "steelblue", size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Densidad de Distribución de Grado log-log",
       x = "Grado (k) - log", y = "Densidad P(k) - log")


##MODELO NULOS 
#Parámetros
n <- vcount(g)
k <- round(mean(degree(g)))
p <- edge_density(g)
n_boot <- 500

get_degrees <- function(g) degree(g)

#Bootstrapping
deg_er_list <- replicate(n_boot, get_degrees(erdos.renyi.game(n, p.or.m = p, type="gnp")), simplify = FALSE)
deg_ba_list <- replicate(n_boot, get_degrees(barabasi.game(n, m=round(k/2), directed=FALSE)), simplify = FALSE)

dens_summary_sd <- function(deg_list){
  all_deg <- unlist(deg_list)
  x_seq <- seq(min(all_deg), max(all_deg), length.out=200)
  
  dens_list <- lapply(deg_list, function(d) {
    d_res <- density(d, from=min(x_seq), to=max(x_seq), n=length(x_seq))
    approx(d_res$x, d_res$y, xout = x_seq)$y
  })
  
  dens_mat <- do.call(rbind, dens_list)
  
  data.frame(
    grado = x_seq,
    mean = apply(dens_mat, 2, mean),
    sd = apply(dens_mat, 2, sd)
  )
}

df_er <- dens_summary_sd(deg_er_list) %>% mutate(red="Erdős–Rényi")
df_ba <- dens_summary_sd(deg_ba_list) %>% mutate(red="Barabási–Albert")

#Obs
deg_obs <- degree(g)
dens_obs <- density(deg_obs, n=200)
df_obs <- data.frame(grado=dens_obs$x, mean=dens_obs$y, red="Observada")


df_plot <- bind_rows(
  df_obs,
  df_er %>% mutate(lower=mean-sd, upper=mean+sd),
  df_ba %>% mutate(lower=mean-sd, upper=mean+sd)
)

ggplot(df_plot, aes(x = grado, color = red, fill = red)) +
  geom_density(alpha = 0.2, adjust = 1.5) +
  theme_minimal() + 
  labs(title = "Comparación de distribuciones de grado",
       x = "Grado", y = "Densidad") +
  scale_color_manual(values = c(
    "Observada" = "steelblue",       
    "Erdős–Rényi" = "orange",        
    "Barabási–Albert" = "purple"     
  )) +
  scale_fill_manual(values = c(
    "Observada" = "steelblue",
    "Erdős–Rényi" = "orange",
    "Barabási–Albert" = "purple"
  ))

#log log
ggplot() +
  geom_ribbon(data = df_plot %>% filter(red!="Observada"),
              aes(x=grado, ymin=lower, ymax=upper, fill=red), alpha=0.2) +
  geom_line(data = df_plot, aes(x=grado, y=mean, color=red), size=0.7) +
  theme_minimal() +
  labs(title="Distribución de grado: Observada vs Modelos nulos",
       x="Grado (log)", y="Pk(log") +
  scale_x_log10(limits = c(100,1500)) + 
  scale_y_log10(limits = c(10^-10, 1))+
  scale_color_manual(values=c("Observada"="steelblue", "Erdős–Rényi"="orange", "Barabási–Albert"="purple")) +
  scale_fill_manual(values=c("Erdős–Rényi"="orange", "Barabási–Albert"="purple"))

#COMPARACION DE METRICAS
obs_density <- edge_density(g)
obs_trans <- transitivity(g, type = "global")
n <- vcount(g)
m <- ecount(g)
k <- round(mean(degree(g)))  
p <- obs_density  

#Erdős–Rényi
g_er <- erdos.renyi.game(n = n, p.or.m = p, type = "gnp")
#Barabási–Albert
g_ba <- barabasi.game(n = n, m = round(k/2), directed = FALSE)

#FUNCION DE METRICAS
metrics <- function(g) {
  c(
    dens = edge_density(g),
    cluster = transitivity(g, type = "global")
  )
}

res <- rbind(
  Observada = metrics(g),
  Erdos_Renyi = metrics(g_er),
  Barabasi_Albert = metrics(g_ba)
)
print(res)

#CENSO TRIADICO
p0 <- function(g) sum(degree(g) == 0)

p1 <- function(g) sum(degree(g) == 1)

p2 <- function(g) {
  n <- vcount(g)
  triangles <- p3(g)
  degrees <- degree(g)
  total_edges <- ecount(g)
  p2_count <- sum(choose(degrees,2)) - 3*triangles
  return(p2_count)
}
p3 <- function(g) sum(count_triangles(g)) / 3 

p0_er <- numeric(n_boot)
p1_er <- numeric(n_boot)
p2_er <- numeric(n_boot)
p3_er <- numeric(n_boot)

for(i in 1:n_boot){
  g_sim <- erdos.renyi.game(n, p.or.m = p, type="gnp")
  p0_er[i] <- p0(g_sim)
  p1_er[i] <- p1(g_sim)
  p2_er[i] <- p2(g_sim)
  p3_er[i] <- p3(g_sim)
}

p0_ba <- numeric(n_boot); p1_ba <- numeric(n_boot); p2_ba <- numeric(n_boot); p3_ba <- numeric(n_boot)

for(i in 1:n_boot){
  g_sim <- barabasi.game(n, m=round(k/2), directed=FALSE)
  p0_ba[i] <- p0(g_sim)
  p1_ba[i] <- p1(g_sim)
  p2_ba[i] <- p2(g_sim)
  p3_ba[i] <- p3(g_sim)
}

# Red observada
p0_obs <- p0(g)
p1_obs <- p1(g)
p2_obs <- p2(g)
p3_obs <- p3(g)

df_triads <- data.frame(
  modelo = rep(c("Observada","Erdős–Rényi","Barabási–Albert"), each=4),
  tipo = rep(c("p0","p1","p2","p3"), times=3),
  valor = c(p0_obs, p1_obs, p2_obs, p3_obs,
            mean(p0_er), mean(p1_er), mean(p2_er), mean(p3_er),
            mean(p0_ba), mean(p1_ba), mean(p2_ba), mean(p3_ba))
)

df_triads_filtered <- df_triads %>% filter(tipo %in% c("p2","p3"))

ggplot(df_triads_filtered, aes(x=tipo, y=valor, fill=modelo)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal() +
  labs(title="Comparación de triadas (p2 y p3)",
       x="Tipo de triada", y="Número de triadas") +
  scale_fill_viridis_d()

#ASORTITIVIDAD DE GRADO
deg <- degree(g)
knn_vals <- knn(g)$knn

#Scatterplot
par(mfrow = c(1,1))
plot(deg, knn_vals,
     xlab = "Grado del nodo",
     ylab = "Grado promedio de los vecinos",
     main = "Grado vs. Grado promedio de los vecinos",
     pch = 21, col = "black", bg = "lightblue")

assortativity_degree(g, directed = FALSE) 
##> 0 homofilia (asortativa);< 0 heterofilia (disasortativa).
