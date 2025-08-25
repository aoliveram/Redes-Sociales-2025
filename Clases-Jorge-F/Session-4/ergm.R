rm(list=ls())
require(statnet)  # Otro paquete tipo Igraph. Nota: es conveniente no usar igraph al mismo tiempo por alcance de nombres de funciones
library(netrankr) # paquete para análisis de centralidad
library(intergraph) # convierte grafos entre distintos formatos
library(ergm)
library(texreg) # para generar tablas comparativas de los modelos

vignette("ergm", package="ergm")

############################################
# ERGM 
# Vamos a reproducir el ejemplo visto en clases
# Ejemplo Familias Florencia
############################################

data(florentine_m)
flomarriage <- intergraph::asNetwork(florentine_m)


flomarriage
gplot(flomarriage, main="Matrimonios")
coordinadas_plot <- (gplot(flomarriage))

gplot(flomarriage, main="Matrimonios", coord = coordinadas_plot)

tabla_atributos <- matrix(NA,nrow=16,ncol=2)
for(i in 1:16){
  tabla_atributos[i,1] <- flomarriage$val[[i]]$vertex.names
  tabla_atributos[i,2] <- flomarriage$val[[i]]$wealth
}

############################################
# partamos con el MODELO BASICO (ER) que solo contiene edges
medici1 <- ergm(flomarriage ~ edges)
summary(medici1)
medici1$coefficients

# son n = 16 nodos, por lo tanto: 
16*15/2 # = 120 posibles nodos
(20)/120 # 0.166 es la probabilidad de un nodo
exp(medici1$coefficients)/(1+exp(medici1$coefficients))

# Nótese que esto es equivalente a una regresión logística
y <- sort(as.vector(as.matrix(flomarriage)))[-c(1:16)]
glm(y~1, family=binomial("logit"))

# otra forma de verlo
prob <- mean(y)
log(prob) - log(1-prob)



###########################################
# Consideremos ahora medidas de CLUSTERING 
# La mas basica son las triadas 
medici2 <- ergm(flomarriage ~ edges + triangle)
summary(medici2)
medici2$coef

# para interpretar: log odds 
triadas0 <- medici2$coef[1] # 0 triada
triadas1 <- medici2$coef[1] + medici2$coef[2] # 1 triada
triadas2 <- medici2$coef[1] + 2*medici2$coef[2] # 2 triada

exp(triadas0)/(1+exp(triadas0)) # prob de link con 0 triada
exp(triadas1)/(1+exp(triadas1)) # prob de link con 1 triada
exp(triadas2)/(1+exp(triadas2)) # prob de link con 2 triada

# por lo tanto, queda claro que la probabilidad de conexiones 
# entre dos nodos aumenta con el 
# numero de triada en las que se participa (... pero estadísticamente no mejora el modelo)

############################################
# agreguemos ahora medidas con ATRIBUTOS de los nodos
# en el caso de esta red un atributo importante es riqueza 
summary(flomarriage) # para ver los atributos de los nodos

riqueza <- flomarriage %v% 'wealth' # asigna el atributo wealth a la variable riqueza
gplot(flomarriage, 
      vertex.cex = riqueza/30, 
      main = "Matrimonios segun riqueza \n Red Renacimiento - Florencia",
      coord = coordinadas_plot)

# para agregar un atributo de la red en ergm usamos nodecov('atributo')

medici3 <- ergm(flomarriage ~ edges + triangle + nodecov('wealth'))
summary(medici3)
medici3$coef  # como es de esperar, coeficiente de riqueza es positivo

# en este caso, cada nodo tiene su propio atributo de riqueza por lo que 
# debe considerarse a ambos. 

tabla_atributos
as.sociomatrix(flomarriage)

# veamos la probabilidad que Medicci este conectado con Ridolfi (que lo esta)
# versus que Barbadori este conectado con Ridolfi (que no lo esta) - ignoremos triadas

e_ridolfi_medici <- medici3$coef[1] + medici3$coef[3]*flomarriage$val[[9]]$wealth + 
  medici3$coef[3]*flomarriage$val[[13]]$wealth
e_ridolfi_bardabori <- medici3$coef[1] + medici3$coef[3]*flomarriage$val[[3]]$wealth + 
  medici3$coef[3]*flomarriage$val[[13]]$wealth

exp(e_ridolfi_medici)/(1+exp(e_ridolfi_medici))
exp(e_ridolfi_bardabori)/(1+exp(e_ridolfi_bardabori))


screenreg(list(medici1, medici2, medici3))

############################################
# QUE MAS SE PUEDE HACER
# el paquete ergm trae otras opciones por defecto
# por ejemplo se puede estudiar homofilia con la funcion nodematch()
# para una lista completa de los atributos incluidos
# ver:
search.ergmTerms()
############################################

############################################
# GOOD OF FITNESS
# una vez hecho un modelo, al igual que en un ejercicio de 
# econometria, proyectamos la estimacion para ver que
# tan bien o mal se ajusta a los datos, eso lo hace la funcion
# 'gof()' 
############################################

medici3.gof <- gof(medici3 ~ model) # usando las variables del modelo
medici3.gof
plot(medici3.gof)

# contrastando con otras variables de la red que no estan en el modelo
medici3.gof.degre <- gof(medici3 ~ degree + esp + distance) # esp = edgewise share partners
medici3.gof.degre
plot(medici3.gof.degre)

############################################
# SIMULACIONES DE REDES
# Luego de revisado que el modelo ajusta bien, se crea redes ficticias
# a partir del modelo que deberian ser parecidas a la real
############################################

medici3.sim <- simulate(medici3, nsim = 15)  # nota: 15 simulaciones es poco 
summary(medici3.sim)
par(mfrow=c(1,2)) # ilustramos dos de las simulaciones generadas
gplot(medici3.sim[[1]], coor = coordinadas_plot)
gplot(medici3.sim[[5]], coor = coordinadas_plot)
dev.off()

############################################
# DIAGNOSTICO DE LA ESTIMACION
# Diagnosticar si MCMC fallo o no (degenerative models)
############################################

mcmc.diagnostics(medici3)  # Si no parece aleatorio estamos en problemas


