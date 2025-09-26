## ============================================================
## Difusión en red WS (misma red para ambas simulaciones)
##   1) Determinista (Valente) con rdiffnet(seed.graph = A)
##   2) Kobayashi (incertidumbre) usando exposure() + regla bayesiana
## ============================================================

library(netdiffuseR)
library(igraph)
library(Matrix)

set.seed(1234)

## -------------------------------------------
## Parámetros globales
## -------------------------------------------
n        <- 250            # nodos
k_ws     <- 8              # grado en el anillo base (par)
p_rewire <- 0.10           # prob. de rewire WS
Tmax     <- 30             # horizonte temporal (slices)
seed_p   <- 0.05           # 5% semillas
theta    <- runif(n, min = 0.1, max = 0.5)   # U(0.1, 0.5) (mismos umbrales en ambas corridas)

## Parámetros Kobayashi (elige y compara)
LAMBDA   <- 0.20    # tolerancia a la incertidumbre (menor = más conservador)
BETA     <- 20      # incertidumbre de exposición (mayor = más precisa)

## -------------------------------------------
## (A) Grafo Watts–Strogatz (estático) con netdiffuseR
##     y plot rápido con igraph
## -------------------------------------------
A <- rgraph_ws(n = n, k = k_ws, p = p_rewire, undirected = TRUE)  # matriz de adyacencia
g_ws <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)

plot(g_ws,
     vertex.size = 1, vertex.label = NA,
     layout = layout_with_fr,
     main = sprintf("Watts–Strogatz: n=%d, k=%d, p=%.2f", n, k_ws, p_rewire))
plot(g_ws,
     vertex.size = 1, vertex.label = NA,
     layout = layout_in_circle(g_ws),
     main = sprintf("Watts–Strogatz: n=%d, k=%d, p=%.2f", n, k_ws, p_rewire))

## -------------------------------------------
## (B) Simulación determinista con rdiffnet()
##     Usando la MISMA red A como seed.graph
##     - seed.nodes = "central" (5% mayor grado)
##     - threshold.dist = theta (vector fijo, umbral por nodo)
## -------------------------------------------
dn_det <- netdiffuseR::rdiffnet(
  seed.graph     = A,        # <- MISMO grafo WS pasado directamente
  t              = Tmax,     # recicla A a lo largo del tiempo
  seed.nodes     = "central",
  seed.p.adopt   = seed_p,
  threshold.dist = theta,    # vector de umbrales por nodo
  rewire         = FALSE,
  stop.no.diff   = FALSE,
  name           = "Determinista (Valente) / WS"
)

## Curva de adopción determinista
par(mfrow = c(1,2), mar = c(4,4,2,1))
plot_adopters(dn_det, main = "Determinista (Valente)\nplot_adopters()")

## -------------------------------------------
## (C) Simulación Kobayashi (fuera de rdiffnet)
##     - MISMA red A
##     - MISMO theta
##     - Semillas "central" 5% (idéntica selección)
##     - Exposición con exposure() (normalizada)
##     - Regla bayesiana con LAMBDA (tolerancia) y BETA (precisión)
## -------------------------------------------
deg       <- Matrix::rowSums(A)
n_seeds   <- max(1, round(seed_p * n))
seed_ids  <- order(deg, decreasing = TRUE)[seq_len(n_seeds)]

y0 <- integer(n); y0[seed_ids] <- 1

k_vec <- pmax(Matrix::rowSums(A), 1)  # evita div/0
eps   <- 1e-8

Y_kob <- matrix(0L, n, Tmax + 1)
Y_kob[, 1] <- y0

Tgraph <- length(dn_det$graph)  # debe ser igual a Tmax si lo fijaste así en rdiffnet()
stopifnot(Tgraph == Tmax)  # opcional, para asegurarte

for (t in 1:Tmax) {
  y_t <- Y_kob[, t]
  
  # n x Tgraph (repetimos el estado en todas las columnas exigidas por exposure)
  cumadopt_t <- matrix(as.integer(y_t), nrow = n, ncol = Tgraph)
  
  # NO CAMBIAR esta llamada (según tu requerimiento)
  q_all <- as.numeric(netdiffuseR::exposure(
    graph      = dn_det$graph,
    cumadopt   = cumadopt_t,
    normalized = TRUE
  ))
  
  # Extraer SOLO la exposición del slice t (bloque de longitud n)
  idx_ini <- (t - 1) * n + 1
  idx_fin <- t * n
  q_t <- q_all[idx_ini:idx_fin]
  
  ## Señal ruidosa m~ ~ Bin(k, q)
  q_clip <- pmin(pmax(q_t, eps), 1 - eps)
  mtilde <- rbinom(n, size = k_vec, prob = q_clip)
  
  ## Prior no sesgado por período: alpha = q*BETA/(1-q)
  alpha  <- (q_clip * BETA) / (1 - q_clip)
  
  ## Posterior: Beta(alpha + m~, BETA + k - m~)
  post_a <- alpha + mtilde
  post_b <- BETA  + k_vec - mtilde
  
  ## Regla de activación (Kobayashi):
  F_theta   <- pbeta(theta, shape1 = post_a, shape2 = post_b)
  adopt_now <- (F_theta < LAMBDA)
  
  ## Monotonía (como Valente): quien adopta, queda adoptado
  y_next <- as.integer(adopt_now | (y_t == 1))
  Y_kob[, t + 1] <- y_next
  
  if (all(y_next == y_t)) break
}

## Empaquetar resultado Kobayashi en diffnet (misma red A)
toa_kob <- rep(NA_integer_, n)
for (i in 1:n) {
  hit <- which(Y_kob[i, ] == 1L)
  if (length(hit)) toa_kob[i] <- min(hit)  # t = 1 => semillas
}
dn_kob <- as_diffnet(A, toa = toa_kob)

## Curva de adopción Kobayashi
plot_adopters(
  dn_kob,
  main = sprintf("Kobayashi (incertidumbre)\nλ=%.2f, β=%d", LAMBDA, BETA)
)

## -------------------------------------------
## (D) Resumen comparativo
## -------------------------------------------
final_det <- mean(dn_det$cumadopt[, ncol(dn_det$adopt)] == 1)
final_kob <- mean(Y_kob[, ncol(Y_kob)] == 1)

cat("\n--- RESUMEN ---\n")
cat(sprintf("Determinista (rdiffnet, WS):    adopción final = %.3f\n", final_det))
cat(sprintf("Kobayashi (λ=%.2f, β=%d, WS): adopción final = %.3f\n", LAMBDA, BETA, final_kob))