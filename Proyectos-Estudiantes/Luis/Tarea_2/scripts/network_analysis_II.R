# libs ----
pacman::p_load(
    tidyverse,
    ggplot2,
    ggthemes,
    ggokabeito,
    igraph,
    neuprintr,
    future.apply
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

## descriptive statistics ----

### network size ----

total_vertices <- vcount(graph_obj)
total_vertices

total_edges <- ecount(graph_obj)
total_edges

### subnetwork size ----

vertex_counts <- table(V(graph_obj)$subnetwork)
vertex_counts

# this are edges that start and end within a given subnetwork
edge_sources <- V(graph_obj)$subnetwork[head_of(graph_obj, E(graph_obj))]
edge_targets <- V(graph_obj)$subnetwork[tail_of(graph_obj, E(graph_obj))]

internal_edge_mask <- (edge_sources == edge_targets) & !is.na(edge_sources)
internal_edge_labels <- edge_sources[internal_edge_mask]

edge_counts <- table(internal_edge_labels)
edge_counts

### network diameter ----

# will use an approximation to make it faster, this should represent the lower bound
num_samples <- 150
set.seed(42)

sample_vertices <- sample(V(graph_obj), size = num_samples)

distances_from_sample <- distances(
    graph_obj,
    v = sample_vertices,
    mode = "out"
)

estimated_diameter <- max(
    distances_from_sample[is.finite(distances_from_sample)]
)
estimated_diameter

### subnetwork diameter ----

# subnetworks
unique_labels <- unique(V(graph_obj)$subnetwork)[1:3]
unique_labels

estimated_diameters <- sapply(unique_labels, function(label) {
    subnetwork_vertices <- which(V(graph_obj)$subnetwork == label)
    sub_G <- induced_subgraph(graph_obj, vids = subnetwork_vertices)
    set.seed(42)
    sample_vertices <- sample(V(sub_G), size = num_samples)
    dists <- distances(
        sub_G,
        v = sample_vertices,
        mode = "out"
    )
    max(dists[is.finite(dists)])
})
estimated_diameters

### average path length ----

estimated_apl_full_network <- mean(
    distances_from_sample[is.finite(distances_from_sample)]
)
estimated_apl_full_network

# repeated here just to maintain order
estimated_apl <- sapply(unique_labels, function(label) {
    subnetwork_vertices <- which(V(graph_obj)$subnetwork == label)
    sub_G <- induced_subgraph(graph_obj, vids = subnetwork_vertices)
    set.seed(42)
    sample_vertices <- sample(V(sub_G), size = num_samples)
    dists <- distances(
        sub_G,
        v = sample_vertices,
        mode = "out"
    )
    mean(dists[is.finite(dists)])
})
estimated_apl


### network density ----

full_density <- edge_density(graph_obj)

### subnetwork density ----

subnetwork_densities <- sapply(unique_labels, function(label) {
    subnetwork_vertices <- which(V(graph_obj)$subnetwork == label)
    sub_G <- induced_subgraph(graph_obj, vids = subnetwork_vertices)
    edge_density(sub_G)
})
subnetwork_densities

### network degree distribution ----

in_degrees <- degree(graph_obj, mode = "in")
out_degrees <- degree(graph_obj, mode = "out")

in_degree_dist <- table(in_degrees)
out_degree_dist <- table(out_degrees)

summary(in_degrees)
summary(out_degrees)

### subnetwork degree distribution ----

subnetwork_degrees <- lapply(unique_labels, function(label) {
    subnetwork_vertices <- which(V(graph_obj)$subnetwork == label)
    sub_G <- induced_subgraph(graph_obj, vids = subnetwork_vertices)
    list(
        in_degree = summary(degree(sub_G, mode = "in")),
        out_degree = summary(degree(sub_G, mode = "out"))
    )
})

names(subnetwork_degrees) <- unique_labels
subnetwork_degrees

### network clusterization ----

avg_clustering <- transitivity(graph_obj, type = "average")

### subnetwork clusterization ----

subnetwork_clustering <- sapply(unique_labels, function(label) {
    subnetwork_vertices <- which(V(graph_obj)$subnetwork == label)
    sub_G <- induced_subgraph(graph_obj, vids = subnetwork_vertices)
    transitivity(sub_G, type = "average")
})
subnetwork_clustering

### centrality metrics ----

### eigenvector ----

subnetwork_eigen <- lapply(unique_labels, function(label) {
    subnetwork_vertices <- which(V(graph_obj)$subnetwork == label)
    sub_G <- induced_subgraph(graph_obj, vids = subnetwork_vertices)
    eigen_centrality(sub_G, directed = TRUE)$vector
})

names(subnetwork_eigen) <- unique_labels
subnetwork_eigen

#### plot ----
p1_df <- enframe(
    subnetwork_eigen,
    name = "subnetwork",
    value = "centrality"
) %>%
    unnest(centrality) %>%
    filter(centrality != 0)


p1 <- p1_df %>%
    ggplot(aes(
        subnetwork, centrality,
        fill = subnetwork
    )) +
    geom_violin(
        trim = TRUE,
        alpha = 0.7
    ) +
    geom_boxplot(
        width = 0.1,
        fill = "white",
        outlier.shape = NA
    ) +
    scale_y_continuous(
        transform = "log10"
    ) +
    xlab("Brain region") +
    ylab("Eigen vector centrality score (log 10)") +
    ggpubr::theme_pubr() +
    theme(legend.position = "none")
p1

### betweennes centrality ----

# removes node with low synaptic weight to speed things up
threshold <- 5
G_simplified <- subgraph.edges(graph_obj, E(graph_obj)[weight >= threshold])
betweenness_simplified <- betweenness(G_simplified, directed = TRUE)

betweenness_df <- tibble(
    subnetwork = V(G_simplified)$subnetwork,
    betweenness = betweenness_simplified
)

bet_mdl <- lm(
    data = betweenness_df %>%
        filter(
            betweenness != 0,
            subnetwork %in% c("cx", "mb", "sez")
        ),
    log10(betweenness) ~ subnetwork
)
summary(bet_mdl)

emmeans::emmeans(
    bet_mdl,
    revpairwise ~ subnetwork,
    type = "response"
)

p2 <- betweenness_df %>%
    filter(
        betweenness != 0,
        subnetwork %in% c("cx", "mb", "sez")
    ) %>%
    ggplot(aes(
        subnetwork, betweenness,
        fill = subnetwork
    )) +
    geom_violin(
        trim = TRUE,
        alpha = 0.7
    ) +
    geom_boxplot(
        width = 0.1,
        fill = "white",
        outlier.shape = NA
    ) +
    scale_y_continuous(
        transform = "log10"
    ) +
    xlab("Brain region") +
    ylab("Betweenness centrality (log 10)") +
    ggpubr::theme_pubr() +
    theme(legend.position = "none")
p2

## null models -----

### ER ----

n_sim <- 100

plan(multisession)

null_metrics <- future_replicate(n_sim, {
    g_er <- sample_gnm(
        total_vertices,
        total_edges,
        directed = TRUE
    )
    sample_nodes <- sample(V(g_er), size = 150)
    dists <- distances(g_er, v = sample_nodes)
    apl_er_approx <- mean(dists[is.finite(dists)])
    c(
        cc = transitivity(g_er, type = "average"),
        apl = apl_er_approx
    )
})
plan(sequential)

cc_null_dist <- null_metrics["cc", ]
apl_null_dist <- null_metrics["apl", ]

#### statistical comparison -----

p_val_cc <- sum(cc_null_dist >= avg_clustering) / n_sim
p_val_cc

p_val_apl <- sum(apl_null_dist >= estimated_apl_full_network) / n_sim
p_val_apl


### SM ----
total_degrees <- in_degrees + out_degrees
avg_total_degree <- round(mean(total_degrees) / 2)

g_lattice <- sample_smallworld(
    dim = 1,
    size = n_real,
    nei = avg_total_degree,
    p = 0
)

gamma <- avg_clustering / mean(cc_null_dist)
lambda <- estimated_apl_full_network / mean(apl_null_dist)
sigma <- gamma / lambda
sigma

### BA -----
m_param <- round(m_real / n_real)

g_ba <- sample_pa(
    n_real,
    power = 1,
    m = m_param,
    directed = TRUE
)

degree_ba <- degree(g_ba, mode = "in")

df_real <- as_tibble(in_degree_dist)
names(df_real) <- c("degree", "frequency")
df_real$network <- "Real network"

df_ba <- as_tibble(table(degree_ba))
names(df_ba) <- c("degree", "frequency")
df_ba$network <- "BA model"

ba_plot_df <- bind_rows(df_real, df_ba) %>%
    mutate(
        degree = as.numeric(as.character(degree))
    )


p3 <- ba_plot_df %>%
    ggplot(aes(
        degree, frequency,
        color = network
    )) +
    geom_point(alpha = 0.7) +
    scale_x_log10() +
    scale_y_log10() +
    ggpubr::theme_pubr() +
    xlab("In-degree (k)") +
    ylab("Frequency P(k)")
p3

p4 <- df_real %>%
    mutate(degree = as.numeric(degree)) %>%
    ggplot(aes(
        degree, frequency
    )) +
    geom_point(
        alpha = 0.7,
        color = "dodgerblue"
    ) +
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(
        method = "lm",
        se = FALSE,
        color = "red",
        linetype = "dashed"
    ) +
    xlab("In-degree (k)") +
    ylab("Frequency P(k)") +
    ggpubr::theme_pubr()
p4
