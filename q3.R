# Q3: ALGORITHM
#
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(RcppAlgos) # For comboGrid function


set.seed(1234)

# Make the graph from Q2
g1 = c(1,2,8,9,7)
g2 = (1:9)[-g1]
df = data.frame(t(matrix(c(1,2,1,8,1,9,1,7,1,3,1,6,3,4,3,6,4,6,6,5), nrow=2)))
graph = graph_from_data_frame()
# Assign group to each node
V(graph)$group = sapply(V(graph), function(i) ifelse(i %in% g1, 1, 2))


# Plot
ggraph(graph, layout='fr') + 
    geom_edge_link(color='#4f555e') + 
    geom_node_point(aes(color=factor(group)), size=6) +
    theme_void()

random_partition = function(vertices) {
    # Return a list of two vectors that partition the input vector
    # TODO maybe change this to return indices of vectors per partition?

    # Toss a coin for every vertex whether to include them in set A or B.
    x = rbinom(length(vertices), 1 , 1/2)
    partition = list(A=vertices[x==0], B=vertices[x == 1])
    return(partition)
}

edges_between_groups = edges_between_groups(u, v, graph) {
    # Get's E_{uv} term
    V(graph)[group==u || group==v] %>% induced_subgraph()
                                   %>% gsize()
}

loglikelihood = function(graph) {
    m = gsize(graph) # Get number of edges of the graph
    c = V(graph)$group %>% unique() # TODO this is not such a good way of doing this
    grid = comboGrid(1:c, 1:c, repetition=F) # Make a 2 column grid of combinations of groups
    grid$E = map2(grid[,1], grid[,2], ~ edges_between_groups(.x, .y, graph)) # Map over combinations and get E_uv term
    grid$kappa = map2(grid[,1], grid[,2], ~ ...) # TODO 
    sum(
        pmap(list(grid$E, grid$kappa1, grid$kappa2),
             function(E, kappa1, kappa2) (E/m) * log( (E/2*m) * ( (kappa1/2*m)*(kappa2/2*m) ))
        ) %>% unlist()
    )
}

possible_moves = function(G, partition, frozen_indices) {
    # Get's possible moves for a given Graph, unfrozen nodes and partition
    ...
}

makeAMove = function(G, partition, c, frozen_indices) {
    # Would this function be better if it was recursive?
    for (node in V(G)[-frozen_indices]) {
        for (c in possible_moves(G, partition, frozen_indices)) {
                 ...
        }
    }
}
