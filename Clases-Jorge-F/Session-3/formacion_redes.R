rm(list=ls())
require(igraph)
require(dplyr)
require(here)
aqui <- here()

a <- 5 
prob <- 0.10
matriz0 <- rbinom(size = 1, n = a ^ 2,   p = prob) %>% matrix(., nrow = a, ncol = a)  # size=1 implica que la matriz tendrá 1s y 0s
red0 <- graph_from_adjacency_matrix(matriz0)
l.circulo <- layout_in_circle(red0)
plot(red0,edge.arrow.size=0, layout = l.circulo)

a1 <- 10
prob <- 0.15
matriz0 <- rbinom(size = 1, n = a ^ 2,   p = prob) %>% matrix(., nrow = a, ncol = a)  
m <- (a1-1)*prob
red0 <- graph_from_adjacency_matrix(matriz0)
grado <- degree_distribution(red0)
plot(density(grado), main = paste0(a1, " nodos"), xlim = c(0,0.5), ylim = c(0,6))

for(x in 1:40){
  a <- a1 + 10*x
  p2 <- m/a
  matriz0 <- rbinom(size = 1, n = a ^ 2,   p = p2) %>% matrix(., nrow = a, ncol = a)  
  red0 <- graph_from_adjacency_matrix(matriz0)
  grado <- degree_distribution(red0)
  plot(density(grado), main = paste0(a, " nodos"), xlim = c(0,0.5), ylim = c(0,6))
}


###################################################################################
###################################################################################
# 1. Generamos una red aleatoria a partir de una matriz de adyacencia asumiendo
# una distribucion binomial
###################################################################################
###################################################################################
set.seed(1234)

a <- 1000 
prob <- 0.15
matriz <- rbinom(size = 1, n = a ^ 2,   p = prob) %>% matrix(., nrow = a, ncol = a) 
red <- graph_from_adjacency_matrix(matriz)
red <- simplify(red, remove.loops = T) 
l.fr <- layout_with_fr(red)
plot(red,
     edge.arrow.size=.001,
     vertex.label=NA,
     vertex.size=10,
     layout=l.fr)
grado <- degree_distribution(red)
t <- paste("Histograma - red aleatoria - ",a," nodos")
plot(density(grado), main = t, xlim=c(0,.15))
hist(grado)

# Lo mismo, pero usando igraph desde el inicio para crear una red aleatoria
# usamos la funcion erdos.renyi.game para crear una red aleatoria, mas adelante la vemos con mas detalle
red2 <- erdos.renyi.game(a, prob, type = "gnp", loops = F) # excluimos loops para hacerla comparable con la anterior
l.fr2 <- layout_with_fr(red)
plot(red2,
     edge.arrow.size=.001,
     vertex.label=NA,
     vertex.size=10,
     layout=l.fr2)
grado2 <- degree_distribution(red2)

par(mfrow=c(2,1))
hist(grado, col = "red")
hist(grado2, col="blue")
dev.off()

###################################################################################
###################################################################################
## 2. Veamos como se compara una red real con una aleatoria de igual distribucion de grado
###################################################################################
###################################################################################

red.e <- read_graph(paste0(aqui,"/data/red_empresas.igraph"))
red.e


# Nos quedamos primero con el componente gigante
componentes.e <- clusters(red.e)
g.e <- which.max(componentes.e$csize) # identificamos el gigante
red.e <- induced.subgraph(red.e, which(componentes.e$membership == g.e)) # nos quedamos con el componente gigante

# vamos a guardar una funcion para extraer el componente gigante en adelante
# asi en el resto del script si necesitamos extraer un componente gigante simplemente
# usamos la funcion nombre_dado_al_componente <- giant.component(el_nombre_de_la_red)

giant.component <- function(graph) {
  cl <- clusters(graph)
  induced.subgraph(graph, which(cl$membership == which.max(cl$csize)))
}

# comparemos la distribucion de grado con una red aleatoria 
# creamos una red con la misma distribucion que red.e pero aleatoria
red.e.degree <- degree(red.e)
red.e.random <- sample_degseq(out.deg = red.e.degree, in.deg = NULL, method = "vl")
plot(degree_distribution(red.e), col="red", type="p")
lines(degree_distribution(red.e.random), col="blue", type = "l")

par(mfrow=c(1,2))
plot(red.e, 
     edge.arrow.size=0,
     vertex.label=NA, 
     vertex.size=degree(red.e)/2, 
     main ="Real")
plot(red.e.random, 
     vertex.label=NA, 
     edge.arrow.size=0,
     vertex.size=degree(red.e.random)/2, 
     main="Aleatoria")
dev.off()

# distancia media
d.e <- round(mean_distance(red.e),3) # por defecto genera media del componente gigante
d.e.random <- round(mean_distance(red.e.random),3)

# clustering medio
t.e <- round(transitivity(red.e, type="global"),3)
t.e.random <- round(transitivity(red.e.random),3) 

resumen.distancias <- matrix(c(d.e,d.e.random,
                               t.e,t.e.random), byrow=T,
                             nrow = 2, ncol = 2)
row.names(resumen.distancias) <- c("distancia empresas", 
                                   "clustering empresas")
colnames(resumen.distancias) <- c("real", "aleatorio")

resumen.distancias # como puede verse, la red aleatoria clusteriza mucho menos que la red real

# Small World Network
??watts.strogatz.game
redsw <- watts.strogatz.game(dim=1, size=15, nei=3, p=0.1, loops = FALSE, multiple = FALSE)
plot(redsw, layout=layout_in_circle(redsw))
distancia_sw <- round(mean_distance(redsw),3)
clustering_sw <- round(transitivity(redsw, type="global"),3)

# Preferential attachment
??sample_pa
redpa <- sample_pa(15,power=1, m=2, directed=F, algorithm="psumtree") 
plot(redpa, layout=layout_in_circle(redpa))
distancia_pa <- round(mean_distance(redpa),3)
clustering_pa <- round(transitivity(redpa, type="global"),3)

################################################################################
################################################################################
# 3. Comparemos ahora la red real con distintos mecanismos de formacion de red
##################################################################################
# dibujamos la red original 
plot(red.e, 
     main="Interlocking", 
     vertex.color="light blue",
     layout=layout.reingold.tilford(red.e, circular=T),
     edge.color="grey",
     edge.width=E(red.e)$weight/10,
     edge.arrow.size=0.1,
     vertex.size=red.e.degree/2, # ve?mos los nodos seg?n centralidad de grado
     vertex.frame.color="blue", 
     vertex.label=NA)

r1 <- table(red.e.degree)
r1 <- as.data.frame(r1)
r1.grado <- as.numeric(r1$red.e.degree)
r1.frec <- r1$Freq
r1.frec 
r1.l <- loess(r1.frec  ~  r1.grado,data=r1)
plot(r1.grado,r1.frec, xlab="grado", ylab="Frecuencia")
lines(r1.l$fitted)
red.e.distancia <- round(mean_distance(red.e),3)
red.e.distancia
red.e.clustering <- round(transitivity(red.e, type="global"),3)
red.e.clustering

plot(r1.grado,r1.frec, xlab="log(grado)", ylab="log(Frecuencia)", log="xy")
segments(10,20,20,1, col="red")
segments(8,20,23,1, col="blue")
text(16,21,"¿Cuál?", col="green")

######## veamos la version Erdos-Renyi de nuestra red ######## 
# erdos.renyi.game(n, p.or.m, type = c("gnp", "gnm"), directed = FALSE, loops = FALSE, ...)
# n: The number of vertices in the graph.
# p.or.m: Either the probability for drawing an edge between two arbitrary vertices (G(n,p) graph), or the number of edges in the graph (for G(n,m) graphs).
# type: The type of the random graph to create, either gnp (G(n,p) graph) or gnm (G(n,m) graph).
# directed: Logical, whether the graph will be directed, defaults to FALSE.
# loops: Logical, whether to add loop edges, defaults to FALSE.
size <- length(V(red.e))
dens <- edge_density(red.e) # probabilidad de un link
er <- erdos.renyi.game(size, dens) # gnp
er.grado <- degree(er)
plot(er, 
     main="Erdös-Renyi", 
     vertex.color="light blue",
     layout=layout.reingold.tilford(er, circular=T),
     edge.color="grey",
     edge.width=E(er)$weight/10,
     edge.arrow.size=0.1,
     vertex.size=er.grado/2, 
     vertex.frame.color="blue", 
     vertex.label=NA)
er.distancia <- round(mean_distance(er),3)
er.clustering <- round(transitivity(er, type="global"),3)

########  veamos la versi?n Strogatz-Watts de nuestra red ######## 
# sample_smallworld(dim, size, nei, p, loops = FALSE, multiple = FALSE)
# dim: Integer constant, the dimension of the starting lattice.
# size: Integer constant, the size of the lattice along each dimension.
# nei: Integer constant, the neighborhood within which the vertices of the lattice will be connected.
# p: Real constant between zero and one, the rewiring probability.
# loops: Logical scalar, whether loops edges are allowed in the generated graph.
# multiple: Logical scalar, whether multiple edges are allowed int the generated graph.
sm <- watts.strogatz.game(1,size,3,0.1) # dim=1, estoy asumiendo vecindarios de 3 nodos y rewiring de 0.1 (se puede mejorar la precisión)
sm.grado <- degree(sm)
plot(sm, 
     main="Strogatz-Watts", 
     vertex.color="light blue",
     layout=layout.reingold.tilford(sm, circular=T),
     edge.color="grey",
     edge.width=E(sm)$weight/10,
     edge.arrow.size=0.1,
     vertex.size=sm.grado/2, 
     vertex.frame.color="blue", 
     vertex.label=NA)
sm.distancia <- round(mean_distance(sm),3)
sm.clustering <- round(transitivity(sm, type="global"),3)

########  veamos la version Barabasi de nuestra red ######## 
# sample_pa(n, power = 1, m = NULL, out.dist = NULL, out.seq = NULL,
#      out.pref = FALSE, zero.appeal = 1, directed = TRUE,
#      algorithm = c("psumtree", "psumtree-multiple", "bag"),
#      start.graph = NULL)
# n; Number of vertices.
# power: The power of the preferential attachment, the default is one, ie. linear
# m: Numeric constant, the number of edges to add in each time step.
# out.dist: Numeric vector, the distribution of the number of edges to add in each time step. This argument is only used if the out.seq argument is omitted or NULL.
# out.seq: Numeric vector giving the number of edges to add in each time step. Its first element is ignored as no edges are added in the first time step.
# out.pref: Logical, if true the total degree is used for calculating the citation probability, otherwise the in-degree is used.
# zero.appeal: The 'attractiveness' of the vertices with no adjacent edges. See details below.
# directed: Whether to create a directed graph.
# algorithm: The algorithm to use for the graph generation.
# start.graph: ... If a graph, then the supplied graph is used as a starting graph for the preferential attachment algorithm. 
red.pa <- sample_pa(size,power=1, m=2, directed=F, algorithm="psumtree") 
degree.red <- degree(red.pa)
l <- layout.reingold.tilford(red.pa, circular=T)
plot(red.pa, 
     main="Red - Pref. Attachment", 
     vertex.color="light blue",
     layout=layout.reingold.tilford(red.pa, circular=T),
     edge.color="grey",
     edge.width=E(red.pa)$weight/10,
     edge.arrow.size=0.1,
     vertex.size=degree.red/2, 
     vertex.frame.color="blue", 
     vertex.label=NA)
red.pa.distancia <- round(mean_distance(red.pa),3)
red.pa.clustering <- round(transitivity(red.pa, type="global"),3)


resumen <- matrix(c(red.e.distancia,red.e.clustering, 
  er.distancia, er.clustering,
  sm.distancia,sm.clustering,
  red.pa.distancia,red.pa.clustering),
  nrow=4,ncol=2,byrow=F)
colnames(resumen) <- c("distancia","clustering")
row.names(resumen) <- c("Real","ER","SW","PA")
resumen

################################################################################
################################################################################
#### Un ejercicio mas completo seria hacer varias redes para tener variacion
################################################################################
################################################################################

v <- V(red.e) # listado de vertices
e <- E(red.e) # listado de vinculos

# simular redes con misma cantidad de nodos y densidad
# que la red que estamos analizando

size <- length(v)
total_redes <- 5 # en un ejercicio completo deberían ser +10000 redes
dens <- edge_density(red.e, loops=F)

# matrix para guardar datos, ejemplo: distribucion de grados
resultados <- matrix(NA,nrow=size, ncol=total_redes+1) 
resultados <- as.data.frame(resultados)
colnames(resultados)[1] <- "interlocking"
resultados[,1] <- degree(red.e)


# Dados los resultados, hagámoslo con Small World
for (d in 1:total_redes){
  red <- watts.strogatz.game(1,size,3,dens, loops=FALSE) 
  red <- giant.component(red)
  degree.red <- degree(red)
  l <- layout.reingold.tilford(red, circular=T)
  resultados[,d+1] <- degree(red)
  colnames(resultados)[d+1] <- paste("SW",d,sep="")
}

View(resultados)

plot(density(resultados[,1]), col="red")
lines(density(resultados[,2]), col="skyblue")
lines(density(resultados[,3]), col="skyblue1")
lines(density(resultados[,4]), col="skyblue2")
lines(density(resultados[,5]), col="skyblue3")
lines(density(resultados[,6]), col="skyblue4")

for (d in 1:total_redes){
  red <- sample_pa(n=size,power=1,algorithm = "psumtree") 
  red <- giant.component(red)
  degree.red <- degree(red)
  l <- layout.reingold.tilford(red, circular=T)
  resultados[,d+1] <- degree(red)
  colnames(resultados)[d+1] <- paste("SW",d,sep="")
}

plot(density(resultados[,1]), col="red")
lines(density(resultados[,2]), col="skyblue")
lines(density(resultados[,3]), col="skyblue1")
lines(density(resultados[,4]), col="skyblue2")
lines(density(resultados[,5]), col="skyblue3")
lines(density(resultados[,6]), col="skyblue4")

##################################################################################
## Cómo saber si una red es powerlaw o no? (o fallar en el intento)
##################################################################################

library(poweRlaw)
# recordemos: 
plot(r1.grado,r1.frec, xlab="log(grado)", ylab="log(Frecuencia)", log="xy")
segments(10,20,20,1, col="red")
segments(8,20,23,1, col="blue")
text(16,21,"¿Cuál?", col="green")

# tres candidatos: power law, log normal y poisson
red.e_pl <- displ$new(degree(red.e)) # discrete power law - displ - ajustamos una distribucion power law a los datos
est.e_pl <- estimate_xmin(red.e_pl) # estimamos el valor inferior de la distribucion
red.e_pl$setXmin(est.e_pl) # actualizamos el objeto red.e_pl con la estimacion del minimo

red.e_ln <- dislnorm$new(degree(red.e)) # ajustamos una distribucion log normal - dislnorm 
est.e_ln <- estimate_xmin(red.e_ln) 
red.e_ln$setXmin(est.e_ln)

red.e_po <- dispois$new(degree(red.e)) # ajustamos una distribucion poisson - dispois 
est.e_po <- estimate_xmin(red.e_po) 
red.e_po$setXmin(est.e_po)

plot(red.e_pl, xlim=c(1,40))
lines(red.e_pl, col="red", lwd=2)
lines(red.e_ln, col="green", lwd=2)
lines(red.e_po, col="blue", lwd=2)

# supongamos que por inspeccion visual descartamos poisson, pero no powerlaw ni lognormal
# como antes, necesitamos muestra, no basta un dato
# hacemos bootstrapping

bs <- bootstrap(red.e_pl, no_of_sims=15, threads=2) # 50 es muy poco en una investigacion
plot(bs, trim=0.1) # trim indica el porcentaje de la muestra que no se despliega

# podemos ahora testear si los datos pueden explicarse por una power-law, Ho: Si, por una PL
bs_p <- bootstrap_p(red.e_pl)
plot(bs_p)
bs_p$p  # si es cercano a cero, podemos estar seguros que no es una power law
# si es menor que 0.1. Clauset et al (2009) sugiere que se puede descartar
# Ahora, podria ser el caso que otra distribucion sea mejor, en este caso, la lognormal.
# para comparar ambas, fijamos el mismo punto de partida para ambas
# y calculamos los parametros de la log normal para ese nuevo punto comun

red.e_ln$setXmin(red.e_pl$getXmin())
est.ln <- estimate_pars(red.e_ln)
red.e_ln$setPars(est.ln)

comparacion <-  compare_distributions(red.e_pl, red.e_ln) # Ho: son indistinguibles
comparacion$test_statistic 
# si es positivo, el primer modelo de distribución se comporta mejor
# Si el valor absoluto es grande, más claro es que una distribución es mejor
# Pero el umbral típico es 1.96. Si es menor, se considera que no hay diferencias significativas
