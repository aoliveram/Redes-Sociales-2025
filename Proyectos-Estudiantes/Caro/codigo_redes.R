#install.packages("igraph")
library(igraph)
library(tidyverse)
library(readxl)
library(dplyr)

setwd("D:/OneDrive/DCCS/1er año/Redes Sociales")
datos <- read_excel("Derechos_concedidos2.xlsx")
names(datos)
df <- datos %>%
  rename(
    nombre_solicitante = "Nombre Solicitante",
    tipo_derecho = "Tipo Derecho",
    naturaleza_agua = "Naturaleza del Agua",
    caudal_anual_prom =  "Caudal \r\nAnual\r\nProm \r\n"
  )

df_red <- df %>%
 select("nombre_solicitante", "Cuenca") %>%
 filter(!is.na(nombre_solicitante) & !is.na(Cuenca))

#df_red <- df %>%
 # select("nombre_solicitante", "Cuenca") %>%
  #filter(!is.na(nombre_solicitante) & !is.na(Cuenca)) %>%
  #group_by(nombre_solicitante, Cuenca) %>%
  #summarise(num_derechos = n(), .groups = 'drop')

edges <- df_red
g <- graph_from_data_frame(edges, directed = FALSE)
V(g)$type <- V(g)$name %in% df_red$nombre_solicitante

plot(
  g,
  vertex.color = ifelse(V(g)$type, "red", "lightgreen"),
  vertex.label = NA,
  vertex.size = 5,
  edge.color = "gray80"
)

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

#invisible(utils::memory.limit(64000))
proj <- bipartite_projection(g)
g_solicitantes <- proj$proj1
deg <- degree(g_solicitantes)

plot(g_solicitantes, 
     vertex.size= 3 + deg,
     #vertex.label = NA,
     vertex.label.cex = 0.5,
     vertex.label.dist = 1.2,
     edge.color = "black",
     #edge.width =  E(g_solicitantes)$weight / 2),
     edge.width = 0.7
)

g_cuencas <- proj$proj2
plot(g_cuencas, vertex.size=5, vertex.label = NA)


library(ggplot2)
library(ggraph)
ggraph(g_solicitantes, layout = "fr") + 
  geom_edge_link(color = "gray70", linewidth = 0.3) + 
  geom_node_point(size = 5, color = "lightblue") +    
  geom_node_text(aes(label = name), repel = TRUE, size = 3) + 
  theme_void()
