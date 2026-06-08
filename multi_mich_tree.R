param_init <- function(L, n_leaf, d) {
  post_params <- list()
  n_node <- 2 * n_leaf - 1
  if (L > 0) {
    for (l in 1:L) {
      post_params[[l]] <- list(
        pi_bar = rep(1 / n_node, n_node),
        log_pi_bar = rep(0.0, n_node),
        b_bar = matrix(0.0, nrow = n_node, ncol = d),
        QTb_bar = matrix(0.0, nrow = n_node, ncol = d),
        mu_bar = matrix(0.0, nrow = n_leaf, ncol = d),
        QTmu_bar = matrix(0.0, nrow = n_leaf, ncol = d)
      )
    }
  }
  return(post_params)
}

tree_mich_matrix <- function(
    y, edges, 
    L=0, L_max=length(y)-1, L_auto=FALSE, 
    fit_intercept=TRUE, fit_scale=TRUE, standardize=TRUE,
    tol=1e-5, max_iter=1e4, verbose=FALSE,
    restart=FALSE, increment=1, 
    merge_level=0.95, merge_prob=NULL,
    pi_l=NULL, omega_l=1e-3
  ) {
  
  #### set up ####
  if (!L_auto) L_max <- L
  
  y_raw <- y  # store original y
  
  d <- ncol(y) # dimension of observations
  n_leaf <- nrow(y)  # number of leafs
  n_node <- 2 * n_leaf - 1 # number of internal nodes
  adj_mat <- tree_mat(n_leaf, edges) 
  
  # min prob to keep component when restarting
  keep_level <- 0.9
  
  # detection threshold
  delta <- 0.5
  detect <- ceiling(log(n_node)^(1 + delta))
  
  # merge threshold
  if (is.null(merge_prob)) merge_prob <- detect / n_node^2
  merge_counter <- log(n_node) %/% 2
  
  # min number of increments for auto procedure
  n_search <- max(4, ceiling(log(n_node) / ((1 + restart) * increment)))
  
  # standardize data for numerical stability 
  if (standardize) {
    center <- colMeans(y)
    scale_eigen <- eigen(var(y))
    Q_scale <- scale_eigen$vectors
    lambda_scale <- scale_eigen$values
    if (any(lambda_scale <= 0)) {
      warning("Var(y) is singular. Consider removing collinear columns.")
      lambda_scale <- lambda_scale + 1e-5
    }
    scale <- Q_scale %*% diag(sqrt(lambda_scale)) %*% t(Q_scale)
    inv_scale <- Q_scale %*% diag(1 / sqrt(lambda_scale)) %*% t(Q_scale)
    y <- (y - matrix(center, nrow = n_leaf, ncol = d, byrow = TRUE)) %*% inv_scale
  }
  
  # initialize mu_0
  mu_0 <- rep(0.0, d)
  
  # estimate precision matrix
  y_diff <- diff(y) # difference out mean changes
  
  # remove outliers due to big mean changes
  y_diff_norm <- sqrt(rowSums(y_diff^2)) 
  y_diff <- y_diff[y_diff_norm <= stats::quantile(y_diff_norm, p = 0.75) +  1.5 * stats::IQR(y_diff_norm), ]
  
  # estimate variance
  Sigma <- (t(y_diff) %*% y_diff) / (2 * nrow(y_diff))
  Sigma_eigen <- eigen(Sigma)
  if (any(Sigma_eigen$values <= 0)) {
    warning("Var(y) is singular. Consider removing collinear columns.")
    Sigma_eigen$values <- Sigma_eigen$values + 1e-5
  }
  Q <- Sigma_eigen$vectors
  lambda <- 1 / Sigma_eigen$values
  
  # set prior probabilities
  if (is.null(pi_l)) pi_l = matrix(1 / n_node, nrow = n_node, ncol = max(1, L))
  if (is.array(pi_l)) log_pi_l = log(pi_l)
  else if (pi_l == "weighted") {
    log_pi_l = log_tree_mean_prior(d, edges)
    log_pi_l = sapply(1:L, function(l) log_pi_l)
  }
  
  # initializing posterior parameters
  post_params <- param_init(L, n_leaf, d)
  
  #### fit model and merge components ####
  merged <- FALSE # flag indicating components have been merged
  while (!merged) {
    fit <- tree_mich_matrix_cpp(
      y, edges, L, 
      mu_0, lambda, Q, fit_intercept, fit_scale, 
      tol, max_iter,
      verbose = verbose & !L_auto,
      log_pi_l, omega_l,
      post_params
    )
    
    merged <- TRUE

    # identify components to merge
    if (L > 1) {
      # extract posterior change-point probabilities
      pi_bar_l <- sapply(1:L, function(l) post_params[[l]][["pi_bar"]])
      
      # only merge columns with credible sets with length less than detect
      cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
      keep <- sapply(cred_sets, length) <= detect
      keep_mat <- matrix(FALSE, ncol = L, nrow = L)
      keep_mat[keep, keep] <- TRUE
      diag(keep_mat) <- FALSE
      
      # compute pairwise merge probabilities
      merge_prob_mat <- t(pi_bar_l) %*% pi_bar_l
      diag(merge_prob_mat) <- 0
      
      # merge components with pairwise merge probabilities > merge_prob
      while (L > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
        merged <- FALSE
        L <- L - 1
        
        # identify components with largest pairwise merge probabilities
        merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]
        
        # keep probabilities of component with largest posterior probability
        if (max(pi_bar_l[,merge_dex[2]]) > max(pi_bar_l[,merge_dex[1]])) {
          pi_bar_l[,merge_dex[1]] <- pi_bar_l[,merge_dex[2]]
          post_params[[merge_dex[1]]][["pi_bar"]] <- post_params[[merge_dex[2]]][["pi_bar"]]
          post_params[[merge_dex[1]]][["log_pi_bar"]] <- post_params[[merge_dex[2]]][["log_pi_bar"]]
          merge_prob_mat[merge_dex[1], ] <- merge_prob_mat[merge_dex[2], ]
          merge_prob_mat[, merge_dex[1]] <- merge_prob_mat[, merge_dex[2]]
        }
        
        # store merged mean parameters
        post_params[[merge_dex[1]]][["QTb_bar"]] <- post_params[[merge_dex[1]]][["QTb_bar"]] + post_params[[merge_dex[2]]][["QTb_bar"]]
        
        # drop merged component
        post_params[merge_dex[2]] <- NULL
        pi_bar_l <- pi_bar_l[,-merge_dex[2]]
        merge_prob_mat <- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
        keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
        if (L_auto) {
          log_pi_l <- log_pi_l[,-merge_dex[2], drop=FALSE]
        }
      }
    }
    if (verbose & !merged) print(paste0("Merging to L = ", L))
  }
  
  # if components were merged out use auto procedure to increase to desired L
  merge_flag <- (L < L_max & !L_auto)
  merge_elbo <- -Inf
  if (merge_flag) increment <- 1 # ensures we don't overshoot number of components
  
  #### auto procedure ####
  
  # 1. Increase L, if ELBO increases, set restart == TRUE and increase L again
  # 2. If ELBO decreases, decrease counter and fit another component
  # 3. If ELBO decreases, counter == 0, and restart == TRUE, reset components
  #    to null model (except for this with concentrated probs) and refit
  # 4. If ELBO decreases, counter == 0, restart == FALSE, and merge_counter > 0,
  #    set fit to best model so far and merge components if no merges return
  #    model, otherwise decrease merge_count and go back to 1
  # 5. If merge_count = 0 return fit
  
  last_restart <- ifelse(restart, 2, Inf)
  
  if (L_auto | merge_flag) {
    counter <- ifelse(merge_flag, L_max-L, n_search) # number of searches after max elbo
    elbo <- fit$elbo[length(fit$elbo)] # store current value of ELBO
    elbo_new <- elbo
    
    if (verbose) print(paste0("L = ", L, ": ELBO = ", elbo_new))
    
    # continue search until n_search exhausted or max components exceeded
    while (L < L_max) {
      # increment dimension of parameters
      for (i in 1:increment) {
        post_params[[L+i]] <- param_init(1, n_leaf, d)[[1]]
      }
      if (L > 0) {
        post_params <- post_params[c((L+1):(L+increment), 1:L)]
        if (L_auto) log_pi_l <- cbind(matrix(log_pi_l[,1], nrow = n_node, ncol = increment), log_pi_l)
        else log_pi_l <- log_pi_l[, c((L_max - L + 1):L_max, 1:(L_max - L))]
      }
      L <- L + increment
      
      # fit incremented model
      fit_new <- tree_mich_matrix_cpp(
        y, edges, L, 
        mu_0, lambda, Q, fit_intercept, fit_scale, 
        tol, max_iter,
        verbose = FALSE,
        log_pi_l, omega_l,
        post_params
      )
      
      # test if model improved or merge/restart ####
      elbo_new <- fit_new$elbo[length(fit_new$elbo)]
      if (verbose) print(paste0("L = ", L, ": ELBO = ", elbo_new,"; Counter: ", counter))
      
      if (elbo_new > elbo | merge_flag) {
        elbo <- elbo_new
        fit <- fit_new
        counter <- ifelse(merge_flag, L_max-L, n_search)
        if (last_restart < Inf) restart <- TRUE
      } else {
        counter <- counter - 1
      }
      
      if (counter == 0 | L == L_max) {
        if (restart) {
          if (verbose) print(paste0("Restarting at L = ", L))
          restart <- FALSE
          last_restart <- L
          
          pi_bar_l <- sapply(1:L, function(l) post_params[[l]][["pi_bar"]])
          
          # reorder by longest blocks first
          chp <- apply(pi_bar_l, 2, which.max)
          block_order <- order(colSums(adj_mat[,chp]))
          post_params <- post_params[block_order]
          pi_bar_l <- sapply(1:L, function(l) post_params[[l]][["pi_bar"]])
          
          # identify components with max prob > keep_level
          keep <- apply(pi_bar_l, 2, max) > keep_level
          
          # reset components with max prob < keep_level to null components
          if (sum(keep) < L) {
            keep_inc <- max(sum(!keep) - increment, 0)
            L <- sum(keep) + keep_inc
            L <- L - increment
            
            log_pi_l <- sapply(1:L, function(i) log_pi_l[,1, drop = FALSE])
            post_params <- post_params[c(which(!keep)[seq_len(keep_inc)], which(keep))]
            for (l in seq_len(keep_inc)) {
              post_params[[l]] <- param_init(1, n_leaf, d)[[1]]
            }
          }
          counter <- ifelse(merge_flag, L_max-L, n_search)
        } else {
          # merging components
          merge_counter <- merge_counter - 1
          no_merges <- TRUE # flag indicating no components were merged
          merged <- FALSE  # flag indicating components have been merged
          if (verbose) print(paste0("Merging. Merge Counter: ", merge_counter))
          
          # set fit to best model so far
          L <- fit$L
          post_params <- fit$post_params
          
          while (!merged) {
            fit <- tree_mich_matrix_cpp(
              y, edges, L, 
              mu_0, lambda, Q, fit_intercept, fit_scale, 
              tol, max_iter,
              verbose = FALSE,
              log_pi_l, omega_l,
              post_params
            )
            
            merged <- TRUE
            
            # identify components to merge
            if (L > 1) {
              # extract posterior change-point probabilities
              pi_bar_l <- sapply(1:L, function(l) post_params[[l]][["pi_bar"]])
              
              # only merge columns with credible sets with length less than detect
              cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
              keep <- sapply(cred_sets, length) <= detect
              keep_mat <- matrix(FALSE, ncol = L, nrow = L)
              keep_mat[keep, keep] <- TRUE
              diag(keep_mat) <- FALSE
              
              # compute pairwise merge probabilities
              merge_prob_mat <- t(pi_bar_l) %*% pi_bar_l
              diag(merge_prob_mat) <- 0
              
              # merge components with pairwise merge probabilities > merge_prob
              while (L > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
                no_merges <- FALSE
                merged <- FALSE
                L <- L - 1
                
                # identify components with largest pairwise merge probabilities
                merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]
                
                # keep probabilities of component with largest posterior probability
                if (max(pi_bar_l[,merge_dex[2]]) > max(pi_bar_l[,merge_dex[1]])) {
                  pi_bar_l[,merge_dex[1]] <- pi_bar_l[,merge_dex[2]]
                  post_params[[merge_dex[1]]][["pi_bar"]] <- post_params[[merge_dex[2]]][["pi_bar"]]
                  post_params[[merge_dex[1]]][["log_pi_bar"]] <- post_params[[merge_dex[2]]][["log_pi_bar"]]
                  merge_prob_mat[merge_dex[1], ] <- merge_prob_mat[merge_dex[2], ]
                  merge_prob_mat[, merge_dex[1]] <- merge_prob_mat[, merge_dex[2]]
                }
                
                # store merged mean parameters
                post_params[[merge_dex[1]]][["QTb_bar"]] <- post_params[[merge_dex[1]]][["QTb_bar"]] + post_params[[merge_dex[2]]][["QTb_bar"]]
                
                # drop merged component
                post_params[merge_dex[2]] <- NULL
                pi_bar_l <- pi_bar_l[,-merge_dex[2]]
                merge_prob_mat <- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
                keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
                if (L_auto) {
                  log_pi_l <- log_pi_l[,-merge_dex[2], drop=FALSE]
                }
              }
            }
            if (verbose & !merged) print(paste0("Merging to L = ", L))
          }
          
          elbo <- max(fit$elbo)
          if (elbo > merge_elbo) {
            merge_elbo <- elbo
            fit_merge <- fit
          }
          
          counter <- ifelse(merge_flag, L_max-L, n_search) # reset counter in case searching continues
          
          # return model if no components to merge or reached max number of merge attempts
          if (no_merges | merge_counter == 0) {
            fit <- fit_merge
            break
          }
        } 
      }
    }
  }
  
  #### return model ####
  class(fit) <- "tree.mich.fit"
  fit$y <- y_raw
  
  # calculate correlated parameters
  for (l in seq_len(fit$L)) {
    fit$post_params[[l]][["b_bar"]] <- fit$post_params[[l]][["QTb_bar"]] %*% t(fit$Q) 
    fit$post_params[[l]][["mu_bar"]] <- fit$post_params[[l]][["QTmu_bar"]] %*% t(fit$Q) 
  }
  fit$Sigma <- fit$Q %*% diag(1 / fit$lambda) %*% t(fit$Q)
  
  # rescale data to original units
  if (standardize) {
    # rescale and center posterior parameters
    fit$mu_0 <- c(fit$mu_0 %*% scale) + center
    fit$mu_bar <- matrix(fit$mu_0, nrow = n_leaf, ncol = d, byrow = TRUE)
    for (l in seq_len(fit$L)) {
      fit$post_params[[l]][["b_bar"]] <- fit$post_params[[l]][["b_bar"]] %*% scale
      fit$post_params[[l]][["mu_bar"]] <- fit$post_params[[l]][["mu_bar"]] %*% scale
    }
    fit$Sigma <- scale %*% fit$Sigma %*% scale
  }
  
  # store probs
  fit$pi_bar <- sapply(seq_len(fit$L), function(l) fit$post_params[[l]][["pi_bar"]])
  
  # calculate mean and residual
  fit$mu_bar <- matrix(fit$mu_0, nrow = n_leaf, ncol = d, byrow = TRUE)
  for (l in seq_len(fit$L)) {
    fit$mu_bar <- fit$mu_bar + fit$post_params[[l]][["mu_bar"]]
  }

  return(fit)
}
  