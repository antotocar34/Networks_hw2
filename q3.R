#
#
# Q3: ALGORITHM
#
# Antoine Carnec
# Danilo Mendez Rubio 
# Leon Karreman
#

# This file contains all the code to answert the questions in Q3
# In particular,it implements the DC-SBM algorithm discussed in the lecture notes.
#
#

library(tidyverse)
library(igraph)
library(igraphdata)
library(ggraph)
library(tidygraph)
library(ggthemes)
library(patchwork)
library(latex2exp)
library(ggtext)
library(kableExtra)

## UTILITY FUNCTIONS
##
# Taken from 
# https://stackoverflow.com/questions/45591286/for-r-markdown-how-do-i-display-a-matrix-from-r-variable
write_matex <- function(x, fname="") {
  begin <- "$$\\begin{bmatrix}"
  end <- "\\end{bmatrix}$$"
  X <-
    apply(x, 1, function(x) {
      paste(
        paste(x, collapse = "&"),
        "\\\\"
      )
    })
  cat(begin, X, end)
}

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

random_partition = function(graph, c) {
    # Toss a c-sided dice for every vertex to determine which partition they belong to
    sample(1:c, gorder(graph), replace=T)
}

ll_label = function(likelihood) str_interp("Log-Likelihood: ${round(likelihood,4)}")
## PLOTTING FUNCTIONS
##
# Make the graph from Q2
make_q2_graph = function() {
    g1 = c(1,2,8,9,7)
    edgelist = c(1,2,1,8,1,9,1,7,1,3,1,6,3,4,3,6,4,6,6,5)
    # edgelist=c(1,2, 1, 8)
    df = data.frame(t(matrix(edgelist, nrow=2)))
    graph = tbl_graph(nodes=data.frame(id=1:9), edges=df)
    # Assign group to each node
    V(graph)$group = map_dbl(V(graph), ~ ifelse(.x %in% g1, 1, 2))
    return(graph)
}



plot_graph = function(graph, partition, layout='fr', message="", label="id") {
    graph %>% 
        activate(nodes) %>% 
        mutate(group = partition) %>%
        ggraph(layout=layout, niter=200) + 
        geom_edge_link(color='#4f555e', edge_width=1.2) + 
        geom_node_point(aes(color=factor(group)), size=13) +
        geom_node_text(aes_string(label=label), size=7, col="black") +
        scale_x_reverse() +
        ggtitle(message) +
        labs(color="groups", title=message) +
        theme_void(base_family="PlayFair") +
        theme(plot.title = element_text(lineheight=1.5, size=20, face="bold.italic"),
              legend.position = "none"
              ) +
        scale_color_tableau("Superfishel Stone")
}

plot_graph_fixed_layout = function(layout, partition, message="", label="id") {
    layout %>% 
        mutate(group = partition) %>%
        ggraph() + 
        geom_edge_link(color='#4f555e', edge_width=1.2) + 
        geom_node_point(aes(color=factor(group)), size=13) +
        geom_node_text(aes_string(label=label), size=7, col="black") +
        scale_x_reverse() +
        ggtitle(message) +
        labs(color="groups", title=message) +
        theme_void() +
        theme(plot.title = element_text(lineheight=1.5, size=20, face="bold.italic"),
              legend.position = "none"
              ) +
        scale_color_tableau("Superfishel Stone")
}

## ALGORITHM IMPLEMENTATION
##
likelihood_calc = function(u, v, E_uv, k_u, k_v, m) {
    # Calculate term inside logliklihood summation
    if (E_uv == 0) return(0) # Avoids having to evalutate log(0)
    (E_uv/(2*m)) * log( (E_uv/(2*m)) / ( (k_u/(2*m)) * (k_v/(2*m)) ) )
}

calculate_edges_between_groups = function(graph, partition, c) {
  # 
  V(graph)$group = partition

  graph = graph %>% 
    activate(edges) %>% 
    mutate(frmgrp = .N()$group[from], togrp = .N()$group[to])
  
  E_tibble = graph %>% activate(edges) %>% 
    as_tibble() %>% 
    group_by(frmgrp, togrp) %>% # Group by group combinations and count edges
    count(name="E") %>% 
    full_join(get_combinations(c), by=c("frmgrp", "togrp")) %>% # Join on combinations in case there are unconnected groups
    replace(is.na(.), 0) %>%
    # Paper states the E_uu is twice the number of edges within community u
    mutate(E = if_else(frmgrp == togrp, 2L*E, E))

  return(E_tibble)
}

loglikelihood = function(graph, partition, c) {
    # Calculates log likelihood
    m = gsize(graph) # Get number of edges of the graph

    E_tibble = calculate_edges_between_groups(graph, partition, c)
    
    E = function(u,v) E_tibble %>% 
                      filter(frmgrp==u, togrp==v) %>% 
                      pull(E)
    
    kappa = map_dbl(1:c, ~ sum(degree(graph)[partition==.x]))

    L = sum(
            map2_dbl(
                 E_tibble$frmgrp, E_tibble$togrp, 
                 function(u, v) likelihood_calc(u,v,
                                                E_uv = E(u,v),
                                                k_u=kappa[u], 
                                                k_v=kappa[v],
                                                m
                                                )
            )
        )
    return(L)
}


makeAMove = function(graph, partition, isFrozen, c) {
    graph_frozen = graph %>% activate(nodes) %>% filter( id %in% which(isFrozen!=1) )

    L_max = -Inf
    best_move = NULL
    # Iterate through the non-frozen nodes
    for (node_index in V(graph_frozen)$id) {
        # Iterate through the (c-1) groups for this node
        for (group_index in (1:c)[-partition[node_index]]) {
            new_partition = partition
            new_partition[node_index] = group_index # Make partition where node is in new group
            L = loglikelihood(graph, new_partition, c) # Calculate log likelihood for this partition
            if (L > L_max) {
                L_max = L
                best_move = list(node_index=node_index, new_group=group_index)
            }
        }
    }
    return( list(likelihood=L_max, best_move=best_move) )
}

runOnePhase = function(graph, initial_partition, c) {
    isFrozen = rep(0, gorder(graph)) # Start with no frozen nodes

    partition_history = tibble(partition=list(initial_partition), 
                               likelihood=loglikelihood(graph, initial_partition, c),
                               i=0
                        )

    partition = initial_partition

    for (i in 1:(gorder(graph))) {
        # Given the frozen nodes, go through possible group changes and return best partition
        best_likelihood_and_move = makeAMove(graph, partition, isFrozen, c)
        best_move = best_likelihood_and_move$best_move

        # move the node to it's new group, creating a new partition
        partition[best_move$node_index] = best_move$new_group
        isFrozen[best_move$node_index] = 1 # Set this node as frozen


        # Add the partition to the history
        partition_history = partition_history %>% add_row(likelihood=best_likelihood_and_move$likelihood,
                                                          partition=list(partition),
                                                          i=i)
    }

    max_move = partition_history %>% slice_max(likelihood) %>% slice(1)
    best_partition = max_move$partition %>% unlist()
    #print(max_move)

    return(list( partition=best_partition, 
                 likelihood=max_move$likelihood, 
                 history=partition_history,
                 h=ifelse(identical(initial_partition, best_partition), 1, 0))
    )
}

fitDCBSM = function(graph, c, T) {
    # Make sure graph has id column
    graph = graph %>% activate(nodes) %>% mutate(id = row_number())

    history = tibble(partition=list(), likelihood=numeric(), i=numeric(), phase=numeric())

    h = 0

    initial_partition = random_partition(graph, c)
    best_likelihood = loglikelihood(graph, initial_partition, c)

    current_partition = initial_partition # Track best partition in this variable.

    for (t in 1:T) {
        out = runOnePhase(graph, current_partition, c) # Should we really start at a random partition?

        if (out$h == 1) {
            print("Halting Criterion was met.") 
            break
        }

        # Append to the history
        history = history %>% bind_rows(out$history %>% mutate(phase = t))

        if (out$likelihood > best_likelihood) {
            best_likelihood = out$likelihood
            current_partition = out$partition
        } else {
            break
        }
    }
    return(list(partition=current_partition, 
                likelihood=best_likelihood, 
                history=history,
                h=out$h
                ))
}

omega_matrix_kappa_vector = function(graph, partition, c) {
    # Produce SBM estimate and kappa matrix, as well as adjacency matrix

    E_tibble = calculate_edges_between_groups(graph, partition, c)

    E_matrix = matrix(NA, nrow=c, ncol=c)

    E_matrix[upper.tri(E_matrix, diag=T)] = E_tibble$E
    # Make the matrix symmetric
    E_matrix[lower.tri(E_matrix)] = t(E_matrix[upper.tri(E_matrix)])

    k = degree(graph)

    kappa = map_dbl(1:c, ~ sum( k[partition == .x]) )

    gamma = map_dbl(1:(gorder(graph)), ~ ( k[.x] )/( kappa[partition[.x]] ) )

    order = gorder(graph)
    A = matrix(NA, ncol=order, nrow=order)
    for (i in 1:order) {
        for (j in 1:order) {
            z_i = partition[i]
            z_j = partition[j]
            A[i, j] = gamma[i]*gamma[j]*E_matrix[z_i,z_j]
        }
    }
    
    return(list( kappa=kappa, E_matrix=E_matrix, A=A ))
}

## Q3 (a)

q2_graph = make_q2_graph()
q2_layout <- create_layout(q2_graph, layout="igraph", algorithm="fr")

c = 3
rand_part = random_partition(q2_graph, c)
rand_likelihood = loglikelihood(q2_graph, rand_part, c)
out = makeAMove(q2_graph, rand_part, c(0,1,0,0,0,0,0,0,0), c)

new_partition = rand_part
new_partition[out$best_move$node_index] = out$best_move$new_group

q3a1 = plot_graph_fixed_layout(q2_layout, rand_part, message=ll_label(rand_likelihood))
q3a2 = plot_graph_fixed_layout(q2_layout, partition=new_partition, message=ll_label(out$likelihood))
ggsave("q3a1.png", plot=q3a1)
ggsave("q3a2.png", plot=q3a2)

## Q3 (b)

random_partition_q3ba = random_partition(q2_graph, 3)
q3b1 = plot_graph_fixed_layout(q2_layout, random_partition_q3ba, message=ll_label(loglikelihood(q2_graph, random_partition_q3ba, c=3)))
out = runOnePhase(q2_graph, rand_part, c=3)
q3b2 = plot_graph_fixed_layout(q2_layout, partition=out$partition, message=ll_label(out$likelihood))
ggsave("q3b1.png", plot=q3b1)
ggsave("q3b2.png", plot=q3b2)

q3b3 = out$history %>% 
       mutate(rownum = nrow(.)) %>% 
       ggplot() +
       geom_line(aes(y=likelihood,x=as.integer(1:rownum)), col="steelblue") +
       geom_point(aes(y=likelihood,x=as.integer(1:rownum)), col="steelblue") +
       xlab("") +
       ylab("Log Likelihood") +
       #ggtitle("Evolution of Log Likelihood over phase" ) +
       theme_bw(base_size = 15)

ggsave("q3b3.png", plot=q3b3)

## Q3 (c)

out_q3c = fitDCBSM(q2_graph, c=3, T=25)

q3c1 = plot_graph_fixed_layout(q2_layout, partition=out_q3c$partition, message=ll_label(out_q3c$likelihood))
ggsave("q3c1.png", plot=q3c1)

history_graph = out_q3c$history

report_omega_kappa = function(graph, partition, c) {
    objs = omega_matrix_kappa_vector(graph, partition, c)
    omega_matrix = objs$E_matrix
    kappa_vector = objs$kappa
    write_matex(omega_matrix)
    write_matex(matrix(kappa_vector, nrow=1, ncol=length(kappa_vector)))
}

q3c_plot = history_graph %>% ggplot() +
    geom_line(aes(y=likelihood, x=i, col=factor(phase))) +
    theme_bw()
  

## Q3 (d)

data(karate)
karate = karate %>% as_tbl_graph()

fit_karate_club = function(iterations, c, T=20) {
    history = tibble(iteration=numeric(), partition=list(), likelihood=numeric(), phase=numeric())
    for (j in 1:iterations) {

      print(str_interp("Iteration ${j}"))
      out = fitDCBSM(karate, 2, T)  
      
      history = history %>% bind_rows(out$history %>% mutate(iteration = j))
    }

    best_iteration = history %>% slice_max(likelihood) %>% slice(1)
    return(list(history=history, 
                partition=unlist(best_iteration$partition), 
                likelihood = best_iteration$likelihood)
    )
}

out = fit_karate_club(1, 2)


q3d2 = out$history %>% ggplot() +
                geom_line(aes(y=likelihood, x=i, col=factor(phase), linetype=factor(iteration))) +
                #geom_point(aes(y=likelihood, x=i, col=factor(phase))) +
                ylab("Log-Likelihood") +
                xlab("Within-Phase Node Iteration") +
                scale_color_discrete( name="Phase No.") +
                guides(linetype="none") +
                theme_bw(base_size=15) +
                theme(legend.title = element_text("Phase Number"))
#ggsave(q3d2, file="q3d2.png")

# Set layout for consistent plotting.
karate_layout = create_layout(karate %>% mutate(id=row_number()), layout="igraph", algorithm="fr")

q3d = plot_graph_fixed_layout(karate_layout, out$partition, message=ll_label(out$likelihood), label="label")
ggsave("q3d.png", plot=q3d)


report_omega_kappa(karate, out$partition, 2)

