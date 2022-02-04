#
# Q3: ALGORITHM
#
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)


get_combinations = function(c) {
    # Creates a dataframe of combinations of groups
    x1 = c()
    x2 = c()
    for (i in 1:c) {
        for (j in i:c) {
            x1= append(x1, i)
            x2 = append(x2, j)
        }
    }

    M =cbind(x1, x2)
    Y = as_tibble(M) %>% rename(frmgrp=x1,togrp=x2)
    return(Y)
}

set.seed(1234)


# Make the graph from Q2
g1 = c(1,2,8,9,7)
g2 = (1:9)[-g1]
edgelist = c(1,2,1,8,1,9,1,7,1,3,1,6,3,4,3,6,4,6,6,5)
# edgelist=c(1,2, 1, 8)
df = data.frame(t(matrix(edgelist, nrow=2)))
graph = tbl_graph(nodes=data.frame(id=1:9), edges=df)
# Assign group to each node
V(graph)$group = map_dbl(V(graph), ~ ifelse(.x %in% g1, 1, 2))
partit = V(graph)$group


# Plot the graph from the question
plot_graph = function(graph, partition, message="") {
    graph %>% activate(nodes) %>% mutate(group = partition) %>%
    ggraph(layout='fr') + geom_edge_link(color='#4f555e') + 
    geom_node_point(aes(color=factor(group)), size=5) +
    geom_node_text(aes(label=id)) +
    scale_x_reverse() +
    ggtitle(message) +
    theme_void()
}

random_partition = function(graph, c) {
    # Toss a c-sided dice for every vertex whether to include them in the partition
    sample(1:c, gorder(graph), replace=T)
}

likelihood_calc = function(u, v, E_uv, k_u, k_v, m) {
    # Calculate term inside logliklihood summation
    if (E_uv == 0) return(0) # Avoids having to evalutate log(0)
    (E_uv/(2*m)) * log( (E_uv/(2*m)) / ( (k_u/(2*m)) * (k_v/(2*m)) ) )
}

loglikelihood = function(graph, partition, c) {
    # Calculates log likelihood
    m = gsize(graph) # Get number of edges of the graph
    V(graph)$group = partition

    graph = graph %>% 
        activate(edges) %>% 
        mutate(frmgrp = .N()$group[from], togrp = .N()$group[to])

    E_tibble = graph %>% activate(edges) %>% 
                         as_tibble() %>% 
                         group_by(frmgrp, togrp) %>% # Group by group combinations and count edges
                         count(name="E") %>% 
                         full_join(get_combinations(c), by=c("frmgrp", "togrp")) %>%
                         replace(is.na(.), 0)
    
    E = function(u,v) E_tibble %>% filter(frmgrp==u, togrp==v) %>% pull(E)

    #grid = comboGrid(1:c, 1:c) %>% as.data.frame() %>% rename(g1=Var1, g2=Var2) # Make a 2 column grid of combinations of groups
    #grid$E = map2_dbl(grid[,1], grid[,2], ~ edges_between_groups(.x, .y, graph)) # Map over combinations and get E_uv term
    degrees = degree(graph)
    L = sum(
            map2_dbl(
                 E_tibble$frmgrp, E_tibble$togrp, 
                 function(u, v) likelihood_calc(u,v,
                                                E_uv = E(u,v),
                                                k_u=(sum(degrees[partition==u])), 
                                                k_v=(sum(degrees[partition==v])),
                                                m
                                                )
            )
        )
    return(L)
}

makeAMove = function(graph, partition, isFrozen, c) {
    graph_frozen = graph %>% activate(nodes) %>% filter( id %in% which(isFrozen!=1) )

    L_max = -Inf# loglikelihood(graph, partition, c) # Initial likelihood
    best_move = NULL
    for (node_index in V(graph_frozen)) {
        for (group_index in (1:c)[-partition[node_index]]) {
            print(str_interp("group index: ${group_index}, ${partition[node_index]}"))
            new_partition = partition
            new_partition[node_index] = group_index # Make partition where node is in new group
            L = loglikelihood(graph, new_partition, c) # Calculate log likelihood for this partition
            if (L > L_max) {
                L_max = L
                best_move = list(node_index=node_index, new_group=group_index)
            }
        }
    }
    return( list(L_max=L_max, best_move=best_move) )
}

runOnePhase = function(graph, initial_partition, c) {
    isFrozen = rep(0, gorder(graph))

    partition_history = tibble(partition=list(), likelihood=numeric())

    partition = initial_partition

    for (i in 1:(gorder(graph))) {
        best_likelihood_and_move = makeAMove(graph, partition, isFrozen, c)
        best_move = best_likelihood_and_move$best_move

        # move the node to it's new group, creating a new partition
        partition[best_move$node_index] = best_move$new_group

        # Add the partition to the history
        partition_history = partition_history %>% add_row(likelihood=best_likelihood_and_move$L_max, 
                                                         partition=list(partition)) 

        isFrozen[node_index] = 1 # Set this node as frozen
    }

    max_move = partition_history %>% slice_max(likelihood)

    return(list(partition=max_move$partition, likelihood=max_move$likelihood))
}

fitDCBSM = function(graph, c, T) {
    isFrozen = rep(0, gorder(graph))

    initial_partition = random_partition(graph, c)
    for (t in 1:T) {
        (best_partition, ) = runOnePhase(graph, initial_partition, )
    }
}

c = 3
rand_part = random_partition(graph, c)
rand_likelihood = loglikelihood(graph, rand_part, c)
out = makeAMove(graph, rand_part, c(0,1,0,0,0,0,0,0,0), c)

new_partition = rand_part
new_partition[out$best_move$node_index] = out$best_move$new_group


q3a1 = plot_graph(graph, rand_part, message=str_interp("Likelihood: ${rand_likelihood}"))
q3a2 = plot_graph(graph, partition=new_partition, message=str_interp("Likelihood: ${out$L_max}"))
ggsave("q3a1.png", plot=q3a1)
ggsave("q3a2.png", plot=q3a2)
