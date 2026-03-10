tree_mat <- function(n_leaf, edges) {
  adj <- matrix(0, nrow = n_leaf, ncol = 2 * n_leaf - 1)
  for (i in 1:n_leaf) {
    j <- i
    while (TRUE) {
      adj[i, j] = 1
      if (j == n_leaf + 1) break
      j <- edges[edges[,2] == j, 1]
    }
  }
  return(adj)
}

tree_sim <- function(n, L, min_space = 10, C = 40, sd = 1, detect_window) {
  tree = ape::rtree(n, rooted = TRUE, br = NULL)
  edges <- tree$edge
  adj <- tree_mat(n, edges)
  n_off <- colSums(adj)
  cand <- (n + 2):(2*n - 1)
  cand <- cand[n_off[cand] >= min_space]
  active_nodes <- c()
  active_color <- rep(0, 2*n - 1)
  leaf_color <- adj %*% active_color
  
  while (length(active_nodes) < L) {
    if (length(cand) == 0) return(tree_sim(n, L, min_space))
    else if (length(cand) == 1) node <- cand
    else node <- sample(cand, 1, prob = n_off[cand] / sum(n_off[cand]))
    
    active_color[node] <- max(leaf_color) + 1
    leaf_color <- adj %*% active_color
    if (min(table(leaf_color)) >= min_space) {
      active_nodes <- c(active_nodes, node)
      parent <- edges[edges[,2] == node, 1]
      sibs <- unlist(edges[edges[,1] == parent, 2])
      cand <- setdiff(cand, sibs)
    } else {
      active_color[node] <- 0
      cand <- setdiff(cand, node)
    }
  }
  
  active_count <- rep(0, L)
  for (i in 1:n) {
    j <- i
    while (j != n + 1) {
      if (j %in% active_nodes) {
        active_count[j == active_nodes] <- active_count[j == active_nodes] + 1
        break
      }
      j <- edges[edges[,2] == j, 1]
    }
  }
  
  jumps <- (1 - 2 * rbinom(L,1,0.5)) * rnorm(L, sqrt(C / active_count), sd = 0.1)
  delta <- rep(0, 2*n - 1)
  delta[n+1] <- runif(1, -2, 2)
  delta[active_nodes] <- jumps
  
  mean_signal <- adj %*% delta
  y <- rnorm(n, mean_signal, sd)
    
  return(list(tree = tree, 
              active_nodes = active_nodes,
              delta = delta,
              mean_signal = mean_signal,
              y = y, 
              adj = adj))
}

tip_color <- function(nodes, adj) {
  n <- nrow(adj)
  active_colors <- rep(0, 2*n - 1)
  for (node in nodes) {
    tip_color <- adj %*% active_colors
    active_colors[node] <- max(tip_color) + 1
  }
  return(as.numeric(as.factor(adj %*% active_colors)))
}

window_fn <- function(chp, adj, wp = 0.5) {
  unlist(
    lapply(chp, function(x) {
      y = adj + adj[,x];
      which(colSums(y == 1) / colSums(y >= 1) < wp)
    })
  )
}

tree_hausdorff <- function(est_cp, true_cp, tree) {
  d_tree <- suppressWarnings(dist.nodes(tree))
  max(sapply(true_cp, function(x) min(d_tree[x,est_cp])))
}

