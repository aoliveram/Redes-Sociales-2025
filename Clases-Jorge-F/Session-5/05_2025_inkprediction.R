# predicci?n de links 
rm(list=ls())
require(igraph)
require(igraphdata)

data(karate)
d <- degree(karate)
e <- as.data.frame(get.edgelist(karate, names=F))

coordenadas <- layout.fruchterman.reingold(karate)
plot(karate, layout=coordenadas)

###############################################################################
###############################################################################
###############################################################################
# algoritmos para detectar comunidades
###############################################################################
###############################################################################
###############################################################################

# ver http://igraph.org/c/doc/igraph-Community.html

community_gn <- edge.betweenness.community(karate)   #Girvan-Newman
community_random <- walktrap.community(karate) 
community_fc <- fastgreedy.community(karate)

plot(karate, 
     layout=coordenadas,
     vertex.color=community_random$membership)


###############################################################################
###############################################################################
###############################################################################
# prediccion de links
###############################################################################
###############################################################################
###############################################################################

# calculamos las similaridades sobre la red original  
jaccard <- similarity(karate, vids = V(karate), mode = "all",
                      loops = FALSE, method = "jaccard")
dice <- similarity(karate, vids = V(karate), mode = "all",
                   loops = FALSE, method = "dice")
invlogweighted <- similarity(karate, vids = V(karate), mode = "all",  # Adamic-Adar
                             loops = FALSE, method = "invlogweighted") 

# creamos sets de entrenamiento, prueba y missing

nodos <- length(V(karate)) 
links_potenciales_num <- nodos*(nodos - 1)/2 # posibles links
red_completa <- erdos.renyi.game(nodos,1) # creamos una red completa de igual numero de nodos para tener todas las combinaciones de links posibles
red_completa.df <- as.data.frame(get.edgelist(red_completa))

set_entrenamiento <- e[sample(row.names(e),round(NROW(e)*0.9,0)),]
et <- row.names(set_entrenamiento)
set_prueba <- e[!(rownames(e) %in% c(et))==T,]
ep <- row.names(set_prueba)
set_missing <- red_completa.df[!(rownames(red_completa.df) %in% c(et,ep))==T,]

karate_reducido <- delete_edges(karate,as.list(set_prueba$V1,set_prueba$V2))
plot(karate_reducido,layout=coordenadas)

# ahora que tenemos una red reducida, podemos comparar las m?tricas de similarity y obtener AUC.

auc <- 0
tot <- 0
for(i in 1:NROW(set_prueba)){
    for(j in 1:NROW(set_missing)){
    a1 <- set_prueba$V1[i]
    a2 <- set_prueba$V2[i]
    b1 <- set_missing$V1[j]
    b2 <- set_missing$V2[j]
    

    if(invlogweighted[a1,a2]>invlogweighted[b1,b2]){auc <- auc+1} # si es mayor
    if(invlogweighted[a1,a2]==invlogweighted[b1,b2]){auc <- auc+0.5} # si es igual
    tot <- tot + 1
  }
}

auc/tot

# al repetir el ejercicio con las otras m?tricas se selecciona la que tiene mejor AUC
# en este caso: invlogweighted

# un ejercicio adicional que podemos hacer es estimar la probabilidad de cada link 
# para eso generamos una bbdd nueva con las variables del modelo predictivo (en este caso dos)
# variable dependiente: una dicot?mica con valor 1 si hay conexi?n y 0 si no hay conexi?n
# variable explicativa: la m?trica de similaridad sacada de la matriz invlogweighted
# con eso podemos hacer un probit simple (n?tese que se pueden agregar otros factores tambi?n al modelo)

link_modelo <- matrix(NA,nrow=NROW(red_completa.df), ncol=2)

for(i in 1:NROW(red_completa.df)){
    link_modelo[i,2] <- invlogweighted[red_completa.df$V1[i],red_completa.df$V2[i]]
    ifelse(are_adjacent(karate,
                        red_completa.df$V1[i],
                        red_completa.df$V2[i])==T,
           link_modelo[i]<-1,link_modelo[i]<-0)
  }

link_modelo <- as.data.frame(link_modelo)
modelo_probit <- glm(V1 ~ V2, family = binomial(link = "probit"), 
                data = link_modelo)
summary(modelo_probit)
resultado <- matrix(NA,nrow=NROW(c), ncol=2)
resultado <- as.data.frame(resultado)
resultado <- predict(modelo_probit, link_modelo, type = "response", se.fit = TRUE)
red_completa.df$prediccion  <- resultado$fit
plot(red_completa.df$prediccion)

