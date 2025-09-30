# libs ----
remotes::install_version("Matrix", version = "1.5-3")
pacman::p_load(
    tidyverse,
    ggplot2,
    ggthemes,
    ggokabeito,
    igraph,
    neuprintr,
    future.apply,
    statnet,
    intergraph,
    texreg,
    btergm
)
setwd(this.path::here())

# load graph object (generated in tarea 1)
graph_obj <- read_rds(file = "../data/graph_obj.rds")

# load neuron map
neuron_map <- read_rds(file = "../data/neuron_map.rds")

# add neuron mapping into graph object
ordered_labels <- neuron_map$brain_area[match(
    V(graph_obj)$name,
    neuron_map$neuron_id
)]

V(graph_obj)$subnetwork <- ordered_labels

# filter low degree nodes
total_degrees <- igraph::degree(graph_obj, mode = "all")
nodes_to_keep <- which(total_degrees >= 5)
graph_obj_core <- induced_subgraph(graph_obj, vids = nodes_to_keep)

# convert igraph into statnet object

# get subgraphs first

subsample_size <- 5000

cx_igraph <- induced_subgraph(
    graph_obj_core,
    vids = which(V(graph_obj_core)$subnetwork == "cx")
)
cx_all_nodes <- which(V(graph_obj_core)$subnetwork == "cx")
cx_sample_nodes <- sample(cx_all_nodes, size = subsample_size)
cx_sub_igraph <- induced_subgraph(graph_obj_core, vids = cx_sample_nodes)

mb_igraph <- induced_subgraph(
    graph_obj_core,
    vids = which(V(graph_obj_core)$subnetwork == "mb")
)
mb_all_nodes <- which(V(graph_obj_core)$subnetwork == "mb")
mb_sample_nodes <- sample(mb_all_nodes, size = subsample_size)
mb_sub_igraph <- induced_subgraph(graph_obj_core, vids = mb_sample_nodes)


sez_igraph <- induced_subgraph(
    graph_obj_core,
    vids = which(V(graph_obj_core)$subnetwork == "sez")
)
sez_all_nodes <- which(V(graph_obj_core)$subnetwork == "sez")
sez_sample_nodes <- sample(sez_all_nodes, size = subsample_size)
sez_sub_igraph <- induced_subgraph(graph_obj_core, vids = sez_sample_nodes)


# actual conversion

cx_net <- asNetwork(cx_sub_igraph)
mb_net <- asNetwork(mb_sub_igraph)
sez_net <- asNetwork(sez_sub_igraph)

cx_igraph <- NULL
mb_igraph <- NULL
sez_igraph <- NULL
graph_obj <- NULL
graph_obj_core <- NULL
neuron_map <- NULL

# ergm fitting
library(parallel)
num_cores <- detectCores() - 1

fit_cx <- ergm(
    cx_net ~ edges +
        gwidegree(decay = 0.25, fixed = TRUE) +
        gwesp(decay = 0.25, fixed = TRUE),
    verbose = TRUE,
    estimate = "MPLE",
    control = control.ergm(
        parallel = num_cores,
        parallel.type = "PSOCK"
    )
)


fit_mb <- ergm(
    mb_net ~ edges +
        gwidegree(decay = 0.25, fixed = TRUE) +
        gwesp(decay = 0.25, fixed = TRUE),
    verbose = TRUE,
    estimate = "MPLE",
    control = control.ergm(
        parallel = num_cores,
        parallel.type = "PSOCK"
    )
)

fit_sez <- ergm(
    sez_net ~ edges +
        gwidegree(decay = 0.25, fixed = TRUE) +
        gwesp(decay = 0.25, fixed = TRUE),
    verbose = TRUE,
    estimate = "MPLE",
    control = control.ergm(
        parallel = num_cores,
        parallel.type = "PSOCK"
    )
)

coef(fit_cx)
summary(fit_cx)$coefficients[, "Std. Error"]

# coefficient comparison
compare_coefficients <- function(coef1, coef2, se1, se2) {
    diff <- coef1 - coef2
    se_diff <- sqrt(se1^2 + se2^2)
    z_score <- diff / se_diff
    p_value <- 2 * pnorm(abs(z_score), lower.tail = FALSE)
    return(list("zscore" = z_score, "pval" = p_value))
}

write_rds(x = fit_cx, file = "fit_cx.rds")
write_rds(x = fit_mb, file = "fit_mb.rds")
write_rds(x = fit_sez, file = "fit_sez.rds")

fit_cx <- read_rds(file = "fit_cx.rds")
fit_mb <- read_rds(file = "fit_mb.rds")
fit_sez <- read_rds(file = "fit_sez.rds")

cx_est <- unname(c(coef(fit_cx), summary(fit_cx)$coefficients[, "Std. Error"]))
cx_est

mb_est <- unname(c(coef(fit_mb), summary(fit_mb)$coefficients[, "Std. Error"]))
mb_est

sez_est <- unname(c(coef(fit_sez), summary(fit_sez)$coefficients[, "Std. Error"]))
sez_est

# gwideg statistical comparison
gwideg_mb_cx <- compare_coefficients(
    mb_est[2], cx_est[2],
    mb_est[5], cx_est[5]
)
gwideg_mb_cx

gwideg_cx_sez <- compare_coefficients(
    cx_est[2], sez_est[2],
    cx_est[5], sez_est[5]
)
gwideg_cx_sez

gwideg_mb_sez <- compare_coefficients(
    mb_est[2], sez_est[2],
    mb_est[5], sez_est[5]
)

gwideg_sez_mb <- compare_coefficients(
    sez_est[2], mb_est[2],
    sez_est[5], mb_est[5]
)
gwideg_sez_mb

# gwesp statistical comparison
gwideg_cx_mb <- compare_coefficients(
    cx_est[3], mb_est[3],
    cx_est[6], mb_est[6]
)
gwideg_cx_mb

gwideg_mb_sez <- compare_coefficients(
    mb_est[3], sez_est[3],
    mb_est[6], sez_est[6]
)
gwideg_mb_sez

# euclidean distances
coef_df <- tibble(
    network = c("CX", "MB", "SEZ"),
    gwidegree = c(-16.17, -13.19, -3.76),
    gwesp = c(2.87, 2.85, 1.71)
)

coef_matrix <- as.matrix(coef_df[, c("gwidegree", "gwesp")])
rownames(coef_matrix) <- coef_df$network

distances <- dist(coef_matrix, method = "euclidean")
distances

dist_matrix <- as.matrix(distances)

dist_df <- as.data.frame(as.table(dist_matrix))
colnames(dist_df) <- c("Red1", "Red2", "Distance")

# heatmap plot of distances
p2 <- ggplot(dist_df, aes(x = Red1, y = Red2, fill = Distance)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Distance, 2)), color = "white", size = 5) +
    scale_fill_viridis_c() +
    theme_minimal(base_size = 14) +
    labs(
        title = "Euclidean distance heatmap",
        x = "",
        y = ""
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2

# multivariate parameter comparison
compare_vectors <- function(coef1, se1, coef2, se2) {
    diff <- coef1 - coef2
    cov_matrix <- diag(se1^2 + se2^2)
    chi2_stat <- t(diff) %*% solve(cov_matrix) %*% diff
    df <- length(coef1)
    p_value <- pchisq(as.numeric(chi2_stat),
        df = df,
        lower.tail = FALSE
    )
    list(
        chi2 = as.numeric(chi2_stat),
        df = df,
        p.value = p_value
    )
}

# coeficient vectors
coef_cx <- cx_est[c(2, 3)]
coef_mb <- mb_est[c(2, 3)]
coef_sez <- sez_est[c(2, 3)]
# standard error vectors
se_cx <- cx_est[c(5, 6)]
se_mb <- mb_est[c(5, 6)]
se_sez <- sez_est[c(5, 6)]

# cx vs sez
test_cx_sez <- compare_vectors(
    coef_cx,
    se_cx,
    coef_sez,
    se_sez
)
test_cx_sez

# cx vs mb
test_cx_mb <- compare_vectors(
    coef_cx,
    se_cx,
    coef_mb,
    se_mb
)
test_cx_mb

# mb vs sez
test_mb_sez <- compare_vectors(
    coef_mb,
    se_mb,
    coef_sez,
    se_sez
)
test_mb_sez


# comparative plot

CX_df <- broom::tidy(fit_cx) %>% mutate(network = "CX")
MB_df <- broom::tidy(fit_mb) %>% mutate(network = "MB")
SEZ_df <- broom::tidy(fit_sez) %>% mutate(network = "SEZ")
plot_df <- bind_rows(CX_df, MB_df, SEZ_df) %>%
    mutate(
        term = paste(rep(c("Edges", "Gwideg", "Gwesp"), 3), network, sep = "-"),
        term_ = rep(c("Edges", "Gwideg", "Gwesp"), 3)
    )
plot_df

p1 <- plot_df %>%
    ggplot(aes(
        term, estimate,
        color = term_
    )) +
    geom_pointrange(aes(
        ymax = estimate + std.error,
        ymin = estimate - std.error
    )) +
    facet_wrap(~term_, scale = "free_x") +
    xlab("") +
    ylab("Estimate") +
    ggpubr::theme_pubr() +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
p1





# comparative table
screenreg(list(
    CX = fit_cx,
    MB = fit_mb,
    SEZ = fit_sez
))
