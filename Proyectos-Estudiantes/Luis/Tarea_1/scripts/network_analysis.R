# libs ----
pacman::p_load(
    tidyverse,
    ggplot2,
    ggthemes,
    ggokabeito,
    igraph,
    neuprintr
)
setwd(this.path::here())

# get data ----

## server connection
neuprint_login(
    server = "https://neuprint.janelia.org",
    token = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6Im5pY29uaWNvbHVhcnRlQGdtYWlsLmNvbSIsImxldmVsIjoibm9hdXRoIiwiaW1hZ2UtdXJsIjoiaHR0cHM6Ly9saDMuZ29vZ2xldXNlcmNvbnRlbnQuY29tL2EvQUNnOG9jTGZIWWtfbndCdjBfSUdwMWR5d1BBRVZ6LV9jRDgyRTIxMVlMZFpDajdzdEkxN1pBPXM5Ni1jP3N6PTUwP3N6PTUwIiwiZXhwIjoxOTM0OTI1NTgwfQ.R399P52aj7_d6lc6bDkNs6WuTiXf8V3PzL1X5vNYUZ0",
    dataset = "hemibrain:v1.2.1"
)

## query data
## im interested in the mushroom body, central complex and subesophageal zone

## mushroom body
mb_neurons <- neuprint_find_neurons(
    input_ROIs = "MB(R)",
    output_ROIs = "MB(R)"
)
mb_rois <- c("MB(R)", "MB(L)")
mb_where <- paste0("n.`", mb_rois, "` = TRUE", collapse = " OR ")
mb_cypher_query <- paste("MATCH (n:Neuron) WHERE", mb_where, "RETURN n.bodyId AS bodyid, n.type AS type, n.instance AS instance, n.status AS status")
mb_query_result <- neuprint_fetch_custom(cypher = mb_cypher_query)
mb_neurons <- unlist(purrr::map(mb_query_result$data, 1))


## central complex
## this is a multi-area roi
central_complex_rois <- c("PB", "FB", "EB", "NO")
cx_where <- paste0("n.`", central_complex_rois, "` = TRUE", collapse = " OR ")
cx_cypher_query <- paste("MATCH (n:Neuron) WHERE", cx_where, "RETURN n.bodyId AS bodyid, n.type AS type, n.instance AS instance, n.status AS status")
cx_query_result <- neuprint_fetch_custom(cypher = cx_cypher_query)
cx_neurons <- unlist(purrr::map(cx_query_result$data, 1))

## subesophageal zone
sez_rois <- c("SAD", "GNG")
sez_where <- paste0("n.`", sez_rois, "` = TRUE", collapse = " OR ")
sez_cypher_query <- paste("MATCH (n:Neuron) WHERE", sez_where, "RETURN n.bodyId AS bodyid, n.type AS type, n.instance AS instance, n.status AS status")
sez_query_result <- neuprint_fetch_custom(cypher = sez_cypher_query)
sez_neurons <- unlist(purrr::map(sez_query_result$data, 1))


# neuron map
neuron_map <- tibble(
    neuron_id = c(mb_neurons, cx_neurons, sez_neurons),
    brain_area = c(
        rep("mb", length(mb_neurons)),
        rep("cx", length(cx_neurons)),
        rep("sez", length(sez_neurons))
    )
)
write_rds(x = neuron_map, file = "../../Tarea_2/data/neuron_map.rds")


## neuron ids
neurons_ids <- c(
    mb_neurons,
    cx_neurons,
    sez_neurons
)

## connectivity matrix

## huge matrix, so better to load the object directly
connectivity_matrix <- readRDS(file = "../data/connectivity_matrix.rds")
connectivity_matrix <- neuprint_get_adjacency_matrix(neurons_ids)

saveRDS(
    object = connectivity_matrix,
    file = "../data/connectivity_matrix.rds",
    compress = TRUE
)

# build graph object ----
graph_obj <- graph_from_adjacency_matrix(
    connectivity_matrix,
    mode = "directed",
    weighted = TRUE
)

saveRDS(object = graph_obj, file = "../../Tarea_2/data/graph_obj.rds")

## area map
area_map <- dplyr::bind_rows(
    tibble(bodyid = mb_neurons, area = "MB"),
    tibble(bodyid = cx_neurons, area = "CX"),
    tibble(bodyid = sez_neurons, area = "SEZ")
) %>%
    group_by(bodyid) %>%
    summarise(
        area = paste(sort(unique(area)), collapse = "+")
    )

## add area to graph
graph_vertices <- tibble(bodyid = as.numeric(V(graph_obj)$name))

## join area map to vertices
graph_vertices <- graph_vertices %>%
    left_join(area_map, by = "bodyid")

## set area attribute
V(graph_obj)$area <- graph_vertices$area
area_colors <- c(
    "MB" = "#E69F00",
    "CX" = "#56B4E9",
    "SEZ" = "#009E73",
    "CX+MB" = "#0072B2",
    "CX+SEZ" = "#D55E00",
    "MB+SEZ" = "#CC79A7",
    "CX+MB+SEZ" = "#F0E442"
)
transparent_colors <- sapply(area_colors, function(color) adjustcolor(color, alpha.f = 0.4))
V(graph_obj)$color <- transparent_colors[V(graph_obj)$area]

# plots ----

## to avoid random unconnected neurons
simplified_graph <- delete_edges(graph_obj, E(graph_obj)[weight < 2])
## remove now isolated neurons
simplified_graph <- delete_vertices(
    simplified_graph,
    V(simplified_graph)[degree(simplified_graph) == 0]
)
comps <- components(simplified_graph)
largest_comp_id <- which.max(comps$csize)
final_graph <- induced_subgraph(
    simplified_graph,
    V(simplified_graph)[comps$membership == largest_comp_id]
)

## compute layout
simplified_layout <- layout_with_fr(
    final_graph,
    grid = "grid",
    niter = 1000
)

png("../figures/full_brain_graph.png", width = 5000, height = 5000, res = 300)
plot(
    final_graph,
    layout = simplified_layout,
    vertex.label = NA,
    vertex.size = 2,
    edge.arrow.size = 0.05,
    main = ""
)
legend(
    "bottomleft",
    legend = names(area_colors),
    fill = area_colors,
    bty = "n",
    cex = 1.2
)
dev.off()
