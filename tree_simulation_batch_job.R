source("tree_sim_functions.R")
source("mich_tree.R")
Rcpp::sourceCpp("mich_tree.cpp")

# get job id and set seed
jobid = as.integer(Sys.getenv("SGE_TASK_ID"))
set.seed(jobid)

# change-point packages
cp_pkgs <- c("treeSeg", "ape")
lapply(cp_pkgs, require, character.only = TRUE)

# significance level
level <- 0.9

# MICH sensitivity parameters
tol <- 1e-5
delta <- 0.5
omega_l <- 1e-3

# simulation settings
settings <- data.frame(n = c(rep(100, 2), rep(500, 3), rep(1000, 3)),
                       L = c(2, 5, rep(c(2, 5, 10), 2)),
                       min_space = c(15, 15, rep(30, 3), rep(50, 3)))

cols <- c("method", "n", "L", "min_space", "time", "n_comp", "L_est", 
          "n_detected", "n_covered", "avg_len", "mean_mse",
          "hausdorff_1", "hausdorff_2")

# initialize results matrix
result <- matrix(ncol = length(cols), nrow = 0)
result <- data.frame(result)
names(result) <- cols

for (j in 1:nrow(settings)) {
  print(j)
  # extract settings
  n <- settings$n[j]; L <- settings$L[j]; min_space <- settings$min_space[j];

  # generate data
  tree_data <- tree_sim(n, L, min_space, C = 200, sd = 1) 
  true_cp <- tree_data$active_nodes
  
  # detection threshold
  max_length <- log(n)^(1+delta)
  
  #### MICH Auto ####
  # fit model and time process
  time <- system.time({
    fit <- mich_tree(tree_data$y, tree_data$tree$edge, L_auto = TRUE, tol = tol, 
                     fit_intercept = TRUE, fit_scale = TRUE, standardize = TRUE,
                     max_iter = Inf, restart = TRUE, omega_l = omega_l)
    })[3]
  
  # calculate mean/precision signal MSEs
  mean_mse <- mean((fit$mu_bar - tree_data$mean_signal)^2)

  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$L > 0) {
    post_sets <- mich_sets(fit$pi_bar_l, max_length, level)
    cs <- lapply(post_sets$sets, function(set) set + n)
    est_cp <- post_sets$cp + n
  }
  
  # estimate number of change-points
  L_est <- length(est_cp)
  
  # detected true change-points
  detected <- true_cp[true_cp %in% window_fn(est_cp, tree_data$adj, wp = 0.5)]
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs)))
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))

  # store results
  result[nrow(result) + 1,"method"] <- "MICH Auto"
  result[nrow(result),-1] <- c(n, L, min_space, time, fit$L, L_est, n_detected,
                               n_covered, avg_len, mean_mse,
                               tree_hausdorff(true_cp, est_cp, tree_data$tree, n), 
                               tree_hausdorff(est_cp, true_cp, tree_data$tree, n))
  
  #### MICH Oracle ####
  # fit model and time process
  time <- system.time({
    fit <- mich_tree(tree_data$y, tree_data$tree$edge, L = L, tol = tol, 
                     fit_intercept = TRUE, fit_scale = TRUE, standardize = TRUE, 
                     max_iter = Inf, restart = TRUE, omega_l = omega_l)
    })[3]
  
  # calculate mean/precision signal MSEs
  mean_mse <- mean((fit$mu_bar - tree_data$mean_signal)^2)
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$L > 0) {
    post_sets <- mich_sets(fit$pi_bar_l, max_length, level)
    cs <- lapply(post_sets$sets, function(set) set + n)
    est_cp <- post_sets$cp + n
  }
  
  # estimate number of change-points
  L_est <- length(est_cp)
  
  # detected true change-points
  detected <- true_cp[true_cp %in% window_fn(est_cp, tree_data$adj, wp = 0.5)]
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs)))
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH Ora"
  result[nrow(result),-1] <- c(n, L, min_space, time, fit$L, L_est, n_detected,
                               n_covered, avg_len, mean_mse,
                               tree_hausdorff(true_cp, est_cp, tree_data$tree, n), 
                               tree_hausdorff(est_cp, true_cp, tree_data$tree, n))
  
  #### treeSeg ####
  # fit model and time process
  time <- system.time({
    fit <- treeSeg(tree_data$y, tree_data$tree, alpha = 1 - level, 
                   fam = "gauss", lengths = "dyadic")
    })[3]
  
  # calculate mean/precision signal MSEs
  mean_mse <- mean((fit$mlP - tree_data$mean_signal)^2)
  
  # estimate number of change-points
  L_est <- fit$numbAN
  
  # determine credible sets and change-points
  est_cp <- fit$mlAN
  
  # detected true change-points
  detected <- true_cp[true_cp %in% window_fn(est_cp, tree_data$adj, wp = 0.5)]
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% fit$confSetAN))
  
  # store results
  result[nrow(result) + 1,"method"] <- "treeSeg"
  result[nrow(result),-1] <- c(n, L, min_space, time, NA, L_est, n_detected,
                               n_covered, NA, mean_mse,
                               tree_hausdorff(true_cp, est_cp, tree_data$tree, n), 
                               tree_hausdorff(est_cp, true_cp, tree_data$tree, n))
}

write.csv(result, paste0("./results/result_",jobid,".csv"), row.names = FALSE)
