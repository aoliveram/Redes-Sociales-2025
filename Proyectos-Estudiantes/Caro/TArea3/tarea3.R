library(dplyr)
library(stringr)
library(readxl)

setwd("D:/OneDrive/DCCS/1er año/Redes Sociales/DerechosConcedidos")
datos <- read_excel("DerechosConcedidos.xlsx")

##SEPARACION BASES DE DATOS
#patron1
patron_juridicas <- "(S\\.A\\.|S\\.A|LTDA|LIMITADA|SPA|SOCIEDAD|EMPRESA|COMPAÑ(Í|I)A|
                      CORP|CORPORACI(Ó|O)N|COOP|MUNICIPALIDAD|UNIVERSIDAD|FUNDACION|ASOCIACION|
                      CLUB|AGRÍCOLA|AGRICOLA|SERVICIO|ASOC\\.|EJERCITO|SOC\\.|MINISTERIO|FISCO|MINERA|
                      CODELCO|FERROCARRIL|S\\.C\\.M\\.|CHILE|DIRECCI(Ó|O)N|INC\\.|CIA\\.|CO\\.|CHILI|S\\.|CAJA)"

#patron2
patron_comunidades <- "(COMUNIDAD|OTROS|SUCESION|COMITE|COMITÉ|SIND\\.|JUNTA|SUC\\.|COM\\.)"

datos <- datos %>%
  mutate(
    solicitante_tipo = case_when(
      str_detect(str_to_upper(`Nombre Solicitante`), patron_juridicas) ~ "Persona Jurídica",
      str_detect(str_to_upper(`Nombre Solicitante`), patron_comunidades) ~ "Comunidad",
      TRUE ~ "Persona Natural"
    )
  )
datos <- datos %>%
  rename(
    nombre_solicitante = "Nombre Solicitante",
    caudal_anual = "Caudal \r\nAnual\r\nProm \r\n"
  )

datos <- datos %>%
  select(nombre_solicitante, Cuenca, caudal_anual, Región, solicitante_tipo) %>%
  filter(!is.na(nombre_solicitante), !is.na(Cuenca)) %>%
  distinct()

datos %>% count(solicitante_tipo)

#base de datos entidades juridicas
datos_1 <- datos %>% 
  filter(solicitante_tipo == "Persona Jurídica")

#base de datos personas naturales
dato_2 <- datos %>% 
  filter(solicitante_tipo == "Persona Natural")

df_red_1 <- datos_1 %>%
  arrange(desc("caudal_anual")) %>%  
  slice_head(n = 10000) 
df_red_1 <- df_red_1 %>%
  select(nombre_solicitante, Cuenca, caudal_anual) %>%
  filter(!is.na(nombre_solicitante), !is.na(Cuenca)) %>%
  distinct()

df_red_2 <- datos_2 %>%
  arrange(desc("caudal_anual")) %>%  
  slice_head(n = 10000) 
df_red_2 <- df_red_2 %>%
  select(nombre_solicitante, Cuenca, caudal_anual) %>%
  filter(!is.na(nombre_solicitante), !is.na(Cuenca)) %>%
  distinct()

f_solic_1 <- factor(df_red_1$nombre_solicitante)
f_cuenc_1 <- factor(df_red_1$Cuenca)
f_solic_2 <- factor(df_red_2$nombre_solicitante)
f_cuenc_2 <- factor(df_red_2$Cuenca)

M_1 <- sparseMatrix(
  i = as.integer(f_solic_1),
  j = as.integer(f_cuenc_1),
  x = 1L,
  dims = c(nlevels(f_solic_1), nlevels(f_cuenc_1)),
  dimnames = list(levels(f_solic_1), levels(f_cuenc_1))
)

M_2 <- sparseMatrix(
  i = as.integer(f_solic_2),
  j = as.integer(f_cuenc_2),
  x = 1L,
  dims = c(nlevels(f_solic_2), nlevels(f_cuenc_2)),
  dimnames = list(levels(f_solic_2), levels(f_cuenc_2))
)

g_1 <- graph_from_incidence_matrix(M_1, weighted = TRUE)
g_2 <- graph_from_incidence_matrix(M_2, weighted = TRUE)

layout_bip_1 <- layout_as_bipartite(g_1, hgap = 3, vgap = 1)
layout_bip_2 <- layout_as_bipartite(g_2, hgap = 3, vgap = 1)

col_solicitantes <- adjustcolor("steelblue", alpha = 0.7)
col_cuencas <- adjustcolor("darkorange2", alpha = 0.8)

grados_1 <- degree(g_1)
grados_2 <- degree(g_2)
tam_vertices_1 <- ifelse(V(g_1)$type, 
                         pmax(0.5, pmin(4, sqrt(grados_1/10))), 
                         pmax(1, pmin(6, sqrt(grados_1/5))))
tam_vertices_2 <- ifelse(V(g_2)$type, 
                         pmax(0.5, pmin(4, sqrt(grados_2/10))), 
                         pmax(1, pmin(6, sqrt(grados_2/5))))

if("caudal_anual" %in% names(df_red_1)) {
  df_edges_1 <- df_red_1 %>%
    mutate(
      from = as.integer(f_solic_1),
      to = as.integer(f_cuenc_1) + nlevels(f_solic_1)
    ) %>%
    select(from, to, caudal_anual)
  edge_list_1 <- as_edgelist(g_1)
  for(i in 1:nrow(edge_list_1)) {
    E(g_1)[i]$weight <- 1  
  }
}

if("caudal_anual" %in% names(df_red_2)) {
  df_edges_2 <- df_red_2 %>%
    mutate(
      from = as.integer(f_solic_2),
      to = as.integer(f_cuenc_2) + nlevels(f_solic_2)
    ) %>%
    select(from, to, caudal_anual)
  edge_list_2 <- as_edgelist(g_2)
  for(i in 1:nrow(edge_list_2)) {
    E(g_2)[i]$weight <- 1  
  }
}

rownames(M_1) <- NULL
M_red_1 <- M_1
rownames(M_2) <- NULL
M_red_2 <- M_2

#PROYECCIONES ENTIDADES
W_s_1 <- tcrossprod(M_red_1)
Matrix::diag(W_s_1) <- 0L
w_threshold_1 <- 2L  
W_s_1@x[W_s_1@x < w_threshold_1] <- 0 
W_s_1 <- drop0(W_s_1) 


#PROYECCIONES PERSONAS NATURALES
W_s_2 <- tcrossprod(M_red_2)
Matrix::diag(W_s_2) <- 0L
w_threshold_2 <- 2L  
W_s_2@x[W_s_2@x < w_threshold_2] <- 0 
W_s_2 <- drop0(W_s_2) 

#GRAFOS A PARTIR DE PROYECCIONES
g_solicitantes_1 <- graph_from_adjacency_matrix(W_s_1, mode = "undirected", weighted = TRUE, diag = FALSE)
V(g_solicitantes_1)$name <- as.character(seq_len(vcount(g_solicitantes_1)))

g_solicitantes_2 <- graph_from_adjacency_matrix(W_s_2, mode = "undirected", weighted = TRUE, diag = FALSE)
V(g_solicitantes_2)$name <- as.character(seq_len(vcount(g_solicitantes_2)))

#COMPONENTE MAS GRANDE
#red entidades
is_connected(g_solicitantes_1) 
comp_1 <- components(g_solicitantes_1)
g_gc_s_1 <- induced_subgraph(g_solicitantes_1, which(comp_1$membership == which.max(comp_1$csize)))

#red personas naturales
is_connected(g_solicitantes_2) 
comp_2 <- components(g_solicitantes_2)
g_gc_s_2 <- induced_subgraph(g_solicitantes_2, which(comp_2$membership == which.max(comp_2$csize)))

g1 <- g_gc_s_1
g2 <- g_gc_s_2
deg_1 <- degree(g1)
deg_2 <- degree(g2)
umbral_1 <- 1000
umbral_2 <- 2600

nodos_filtrados_1 <- V(g1)[deg_1 > 0 & deg_1 <= umbral_1]
g_gc_s_filt_1 <- induced_subgraph(g1, vids = nodos_filtrados_1)
g1 <- g_gc_s_filt_1

nodos_filtrados_2 <- V(g2)[deg_2 > 0 & deg_2 <= umbral_2]
g_gc_s_filt_2 <- induced_subgraph(g2, vids = nodos_filtrados_2)
g2 <- g_gc_s_filt_2
is_connected(g1)

deg_1 <- degree(g1)
deg_2 <- degree(g2)

# Gráfico de densidad
df_deg <- data.frame(
  grado = c(deg_1, deg_2),
  tipo  = c(rep("Entidad jurídica", length(deg_1)),
            rep("Persona natural", length(deg_2)))
)

ggplot(df_deg, aes(x = grado, color = tipo, fill = tipo)) +
  geom_density(alpha = 0.4, adjust = 1.5, size = 0.8) +
  theme_minimal() +
  labs(
    title = "Distribución de grado según tipo de solicitante",
    x = "Grado (k)",
    y = "Densidad (P(k))",
    color = "Tipo de solicitante",
    fill = "Tipo de solicitante"
  ) +
  scale_color_manual(values = c("Entidad jurídica" = "steelblue", 
                                "Persona natural" = "orange")) +
  scale_fill_manual(values = c("Entidad jurídica" = "steelblue", 
                               "Persona natural" = "orange")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top"
  )

d_1 <- density(deg_1)
df_d_1 <- data.frame(x = d_1$x, y = d_1$y, tipo = "Entidad jurídica")
d_2 <- density(deg_2)
df_d_2 <- data.frame(x = d_2$x, y = d_2$y, tipo = "Persona natural")
df_all <- rbind(df_d_1, df_d_2)

ggplot(df_all, aes(x = x, y = y, color = tipo)) +
  geom_line(size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() +
  labs(
    title = "Densidad de Distribución de Grado log-log",
    x = "Grado (k) - log",
    y = "Densidad P(k) - log",
    color = "Tipo de solicitante"
  ) +
  scale_color_manual(values = c("Entidad jurídica" = "firebrick",
                                "Persona natural" = "steelblue")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top"
  )

#Métricas
#Tamaño red
num_nodos_1 <- vcount(g1)
num_aristas_1 <- ecount(g1)
num_nodos_2 <- vcount(g2)
num_aristas_2 <- ecount(g2)
#Diámetro
diametro_1 <- diameter(g1, directed = FALSE, unconnected = FALSE)
diametro_2 <- diameter(g2, directed = FALSE, unconnected = FALSE)
#Densidad ---
densidad_1 <- edge_density(g1)
densidad_2 <- edge_density(g2)
#Clusterización ---
#clustering local
clust_local_1 <- transitivity(g1, type = "local", isolates = "zero")
clust_media_1 <- mean(clust_local_1, na.rm = TRUE)
clust_local_2 <- transitivity(g2, type = "local", isolates = "zero")
clust_media_2 <- mean(clust_local_2, na.rm = TRUE)
#Clusterización global
clust_global_1 <- transitivity(g1, type = "global")
clust_global_2 <- transitivity(g2, type = "global")
#Longitud de camino medio
long_camino_medio_1 <- mean_distance(g1, directed = FALSE, unconnected = FALSE)
long_camino_medio_2 <- mean_distance(g2, directed = FALSE, unconnected = FALSE)
#asortatividad
assortativity_degree(g1, directed = FALSE) 
assortativity_degree(g2, directed = FALSE) 

#COMPARACION CON MODELOS NULOS
n_1 <- vcount(g1)
n_2 <- vcount(g2)
k_1 <- round(mean(degree(g1)))
k_2 <- round(mean(degree(g2)))
p_1 <- edge_density(g1)
p_2 <- edge_density(g2)
n_boot <- 500

get_degrees <- function(g) degree(g)

#Bootstrapping
deg_er_list <- replicate(n_boot_1, get_degrees(erdos.renyi.game(n_1, p.or.m = p_1, type="gnp")), simplify = FALSE)
deg_ba_list <- replicate(n_boot_1, get_degrees(barabasi.game(n_1, m=round(k_1/2), directed=FALSE)), simplify = FALSE)
deg_er_list_2 <- replicate(n_boot_2, get_degrees(erdos.renyi.game(n_2, p.or.m = p_2, type="gnp")), simplify = FALSE)
deg_ba_list_2 <- replicate(n_boot_2, get_degrees(barabasi.game(n_2, m=round(k_2/2), directed=FALSE)), simplify = FALSE)

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

#densidad de grado
df_er_1 <- dens_summary_sd(deg_er_list) %>% mutate(red="Erdős–Rényi")
df_ba_1<- dens_summary_sd(deg_ba_list) %>% mutate(red="Barabási–Albert")

df_er_2 <- dens_summary_sd(deg_er_list_2) %>% mutate(red="Erdős–Rényi")
df_ba_2<- dens_summary_sd(deg_ba_list_2) %>% mutate(red="Barabási–Albert")

#Obs
deg_obs_1 <- degree(g1)
dens_obs_1 <- density(deg_obs_1, n=200)
df_obs_1 <- data.frame(grado=dens_obs_1$x, mean=dens_obs_1$y, red="Observada")

deg_obs_2 <- degree(g2)
dens_obs_2 <- density(deg_obs_2, n=200)
df_obs_2 <- data.frame(grado=dens_obs_2$x, mean=dens_obs_2$y, red="Observada")

df_plot_1 <- bind_rows(
  df_obs_1,
  df_er_1 %>% mutate(lower=mean-sd, upper=mean+sd),
  df_ba_1 %>% mutate(lower=mean-sd, upper=mean+sd)
)

df_plot_2 <- bind_rows(
  df_obs_2,
  df_er_2 %>% mutate(lower=mean-sd, upper=mean+sd),
  df_ba_2 %>% mutate(lower=mean-sd, upper=mean+sd)
)

#red1
df_density <- bind_rows(
  data.frame(grado = deg_obs_1, red = "Observada"),
  data.frame(grado = unlist(deg_er_list), red = "Erdős–Rényi"),
  data.frame(grado = unlist(deg_ba_list), red = "Barabási–Albert")
)

ggplot(df_density, aes(x = grado, color = red, fill = red)) +
  geom_density(alpha = 0.1, adjust = 1.5) +  # alpha controla la transparencia
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

#red2
df_density_2 <- bind_rows(
  data.frame(grado = deg_obs_2, red = "Observada"),
  data.frame(grado = unlist(deg_er_list_2), red = "Erdős–Rényi"),
  data.frame(grado = unlist(deg_ba_list_2), red = "Barabási–Albert")
)

ggplot(df_density_2, aes(x = grado, color = red, fill = red)) +
  geom_density(alpha = 0.1, adjust = 1.5) +  # alpha controla la transparencia
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


#metricas
#funcion de metricas
metrics <- function(g) {
  c(
    dens = edge_density(g),
    cluster = transitivity(g, type = "global")
  )
}

metrics_er_1 <- replicate(
  n_boot, 
  metrics(erdos.renyi.game(n_1, p.or.m = p_1, type = "gnp")),
  simplify = TRUE
)
metrics_ba_1 <- replicate(
  n_boot, 
  metrics(barabasi.game(n_1, m = round(k_1/2), directed = FALSE)),
  simplify = TRUE
)

metrics_er_1 <- t(metrics_er_1)
metrics_ba_1 <- t(metrics_ba_1)
colnames(metrics_er_1) <- c("dens", "cluster")
colnames(metrics_ba_1) <- c("dens", "cluster")

metrics_obs_1 <- metrics(g1)

summary_stats <- function(x) {
  c(
    mean = mean(x, na.rm = TRUE),
    lwr = quantile(x, 0.025, na.rm = TRUE),
    upr = quantile(x, 0.975, na.rm = TRUE)
  )
}


er_summary_1 <- apply(metrics_er_1, 2, summary_stats)
ba_summary_1 <- apply(metrics_ba_1, 2, summary_stats)

results_1 <- data.frame(
  Metrica = c("Densidad", "Clustering"),
  Observada = c(metrics_obs_1["dens"], metrics_obs_1["cluster"]),
  ER_Mean = c(er_summary_1["mean","dens"], er_summary_1["mean","cluster"]),
  ER_Lwr  = c(er_summary_1["lwr.2.5%","dens"],  er_summary_1["lwr.2.5%","cluster"]),
  ER_Upr  = c(er_summary_1["upr.97.5%","dens"], er_summary_1["upr.97.5%","cluster"]),
  BA_Mean = c(ba_summary_1["mean","dens"], ba_summary_1["mean","cluster"]),
  BA_Lwr  = c(ba_summary_1["lwr.2.5%","dens"],  ba_summary_1["lwr.2.5%","cluster"]),
  BA_Upr  = c(ba_summary_1["upr.97.5%","dens"], ba_summary_1["upr.97.5%","cluster"])
)

print(results_1, digits = 3)

#densidades
df_density <- bind_rows(
  data.frame(grado = deg_obs_1, red = "Observada"),
  data.frame(grado = unlist(deg_er_list), red = "Erdős–Rényi"),
  data.frame(grado = unlist(deg_ba_list), red = "Barabási–Albert")
)

ggplot(df_density, aes(x = grado, color = red, fill = red)) +
  geom_density(alpha = 0.1, adjust = 1.5) +  # alpha controla la transparencia
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

###2
metrics_er_2 <- replicate(
  n_boot, 
  metrics(erdos.renyi.game(n_2, p.or.m = p_2, type = "gnp")),
  simplify = TRUE
)

metrics_ba_2 <- replicate(
  n_boot, 
  metrics(barabasi.game(n_2, m = round(k_2/2), directed = FALSE)),
  simplify = TRUE
)

metrics_er_2 <- t(metrics_er_2)
metrics_ba_2 <- t(metrics_ba_2)

colnames(metrics_er_2) <- c("dens", "cluster")
colnames(metrics_ba_2) <- c("dens", "cluster")

metrics_obs_2 <- metrics(g2)

er_summary_2 <- apply(metrics_er_2, 2, summary_stats)
ba_summary_2 <- apply(metrics_ba_2, 2, summary_stats)

results_2 <- data.frame(
  Metrica = c("Densidad", "Clustering"),
  Observada = c(metrics_obs_2["dens"], metrics_obs_2["cluster"]),
  ER_Mean = c(er_summary_2["mean","dens"], er_summary_2["mean","cluster"]),
  ER_Lwr  = c(er_summary_2["lwr.2.5%","dens"],  er_summary_2["lwr.2.5%","cluster"]),
  ER_Upr  = c(er_summary_2["upr.97.5%","dens"], er_summary_2["upr.97.5%","cluster"]),
  BA_Mean = c(ba_summary_2["mean","dens"], ba_summary_2["mean","cluster"]),
  BA_Lwr  = c(ba_summary_2["lwr.2.5%","dens"],  ba_summary_2["lwr.2.5%","cluster"]),
  BA_Upr  = c(ba_summary_2["upr.97.5%","dens"], ba_summary_2["upr.97.5%","cluster"])
)

print(results_2, digits = 3)

#densidades
df_density_2 <- bind_rows(
  data.frame(grado = deg_obs_2, red = "Observada"),
  data.frame(grado = unlist(deg_er_list_2), red = "Erdős–Rényi"),
  data.frame(grado = unlist(deg_ba_list_2), red = "Barabási–Albert")
)

ggplot(df_density_2, aes(x = grado, color = red, fill = red)) +
  geom_density(alpha = 0.1, adjust = 1.5) +  # alpha controla la transparencia
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


#censo triadico
#Bootstrapping
deg_er_list_1 <- replicate(n_boot, get_degrees(erdos.renyi.game(n_1, p.or.m = p_1, type="gnp")), simplify = FALSE)
deg_ba_list_1 <- replicate(n_boot, get_degrees(barabasi.game(n_1, m=round(k_1/2), directed=FALSE)), simplify = FALSE)

#funciones
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


p0_er_1 <- numeric(n_boot)
p1_er_1 <- numeric(n_boot)
p2_er_1 <- numeric(n_boot)
p3_er_1 <- numeric(n_boot)

for(i in 1:n_boot){
  g_sim <- erdos.renyi.game(n_1, p.or.m = p_1, type="gnp")
  p0_er_1[i] <- p0(g_sim)
  p1_er_1[i] <- p1(g_sim)
  p2_er_1[i] <- p2(g_sim)
  p3_er_1[i] <- p3(g_sim)
}

p0_ba_1 <- numeric(n_boot)
p1_ba_1 <- numeric(n_boot)
p2_ba_1 <- numeric(n_boot)
p3_ba_1 <- numeric(n_boot)

for(i in 1:n_boot){
  g_sim <- barabasi.game(n_1, m=round(k_1/2), directed=FALSE)
  p0_ba_1[i] <- p0(g_sim)
  p1_ba_1[i] <- p1(g_sim)
  p2_ba_1[i] <- p2(g_sim)
  p3_ba_1[i] <- p3(g_sim)
}

#red1
p0_obs_1 <- p0(g1)
p1_obs_1 <- p1(g1)
p2_obs_1 <- p2(g1)
p3_obs_1 <- p3(g1)

df_triads_1 <- data.frame(
  modelo = rep(c("Observada","Erdős–Rényi","Barabási–Albert"), each=4),
  tipo = rep(c("p0","p1","p2","p3"), times=3),
  valor = c(p0_obs_1, p1_obs_1, p2_obs_1, p3_obs_1,
            mean(p0_er_1), mean(p1_er_1), mean(p2_er_1), mean(p3_er_1),
            mean(p0_ba_1), mean(p1_ba_1), mean(p2_ba_1), mean(p3_ba_1)),
  sd = c(0, 0, 0, 0,                 
         sd(p0_er_1), sd(p1_er_1), sd(p2_er_1), sd(p3_er_1),
         sd(p0_ba_1), sd(p1_ba_1), sd(p2_ba_1), sd(p3_ba_1))
)
df_triads_filtered <- df_triads_1 %>% filter(tipo %in% c("p2","p3"))

ggplot(df_triads_filtered, aes(x=tipo, y=valor, fill=modelo)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal() +
  labs(title="Comparación de triadas (p2 y p3)",
       x="Tipo de triada", y="Número de triadas") +
  scale_fill_viridis_d() +
  geom_errorbar(aes(ymin = valor - sd, ymax = valor + sd),
                width = 0.2, position = position_dodge(0.9))
#red2
p0_er_2 <- numeric(n_boot)
p1_er_2 <- numeric(n_boot)
p2_er_2 <- numeric(n_boot)
p3_er_2 <- numeric(n_boot)

for(i in 1:n_boot){
  g_sim <- erdos.renyi.game(n_2, p.or.m = p_2, type="gnp")
  p0_er_2[i] <- p0(g_sim)
  p1_er_2[i] <- p1(g_sim)
  p2_er_2[i] <- p2(g_sim)
  p3_er_2[i] <- p3(g_sim)
}

p0_ba_2 <- numeric(n_boot)
p1_ba_2 <- numeric(n_boot)
p2_ba_2 <- numeric(n_boot)
p3_ba_2 <- numeric(n_boot)

for(i in 1:n_boot){
  g_sim <- barabasi.game(n_2, m=round(k_2/2), directed=FALSE)
  p0_ba_2[i] <- p0(g_sim)
  p1_ba_2[i] <- p1(g_sim)
  p2_ba_2[i] <- p2(g_sim)
  p3_ba_2[i] <- p3(g_sim)
}

# Red observada
p0_obs_2 <- p0(g2)
p1_obs_2 <- p1(g2)
p2_obs_2 <- p2(g2)
p3_obs_2 <- p3(g2)

df_triads_2 <- data.frame(
  modelo = rep(c("Observada","Erdős–Rényi","Barabási–Albert"), each=4),
  tipo = rep(c("p0","p1","p2","p3"), times=3),
  valor = c(p0_obs_2, p1_obs_2, p2_obs_2, p3_obs_2,
            mean(p0_er_2), mean(p1_er_2), mean(p2_er_2), mean(p3_er_2),
            mean(p0_ba_2), mean(p1_ba_2), mean(p2_ba_2), mean(p3_ba_2)),
  sd = c(0, 0, 0, 0,                 # Observada (si es un valor único, SD=0)
         sd(p0_er_2), sd(p1_er_2), sd(p2_er_2), sd(p3_er_2),
         sd(p0_ba_2), sd(p1_ba_2), sd(p2_ba_2), sd(p3_ba_2))
)

df_triads_filtered_2 <- df_triads_2 %>% filter(tipo %in% c("p2","p3"))


ggplot(df_triads_filtered_2, aes(x=tipo, y=valor, fill=modelo)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal() +
  labs(title="Comparación de triadas (p2 y p3)",
       x="Tipo de triada", y="Número de triadas") +
  scale_fill_viridis_d() 
+ geom_errorbar(aes(ymin = valor - sd, ymax = valor + sd),
                width = 0.2, position = position_dodge(0.9))

#proporciones entre redes
proporcion_p3_jur <- p3_obs_1/n_1
prop_p3_nat <- p3_obs_2/n_2
pro_p2_jur <- p2_obs_1/n_1
prop_p2_nat <- p2_obs_2/n_2
