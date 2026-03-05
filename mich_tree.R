mich_tree <- function(y, edges, 
                      L=0, L_max=length(y)-1, L_auto=FALSE, 
                      fit_intercept=TRUE, fit_scale=TRUE, standardize=TRUE,
                      tol=1e-5, max_iter=1e4, verbose=FALSE,
                      restart=FALSE, increment=1, 
                      merge_level=0.95, merge_prob=NULL,
                      pi_l=NULL, omega_l=1e-3) {

  if (!L_auto) L_max <- L
  
  # number of leafs
  n_leaf <- length(y)
  
  # number of internal nodes
  n_node <- n_leaf - 1
  
  # min prob to keep component when restarting
  keep_level <- 0.9
  
  # detection threshold
  delta = 0.5
  detect <- ceiling(log(n_node)^(1 + delta))
  
  # merge threshold
  if (is.null(merge_prob)) merge_prob <- detect / n_node^2
  merge_counter = log(n_node) %/% 2

  # min number of increments for auto procedure
  n_search <- max(4, ceiling(log(n_node) / ((1 + restart) * increment)))
  
  #### standardize data ####
  if (standardize) {
    center <- mean(y)
    scale <- sd(y)
    y <- (y - center) / scale
  }
  
  #### initialize mu_0 and lambda_0 ####
  mu_0 <- c(0.0)
  lambda_0 <- c(1.0)
  if (fit_scale) {
    # difference out mean changes
    y_diff <- diff(y)
    # remove outliers due to big mean changes
    y_diff <- y_diff[y_diff <= stats::quantile(y_diff, p = 0.75) + 1.5 * stats::IQR(y_diff) &
                       y_diff >= stats::quantile(y_diff, p = 0.25) - 1.5 * stats::IQR(y_diff)]
    lambda_0[1] <- 2.0 / (var(y_diff))
  }
  
  # set prior probabilities
  if (is.null(pi_l)) pi_l = matrix(1 / n_node, nrow = n_node, ncol = max(1, L))
  log_pi_l = log(pi_l)
  
  #### initializing posterior parameters ####
  b_bar_l <- matrix(0.0, nrow = n_node, ncol = L)
  omega_bar_l <- matrix(1.0, nrow = n_node, ncol = L)
  pi_bar_l <- matrix(1 / n_node, nrow = n_node, ncol = L)
  log_pi_bar_l <- matrix(0.0, nrow = n_node, ncol = L)

  #### fit model and merge components ####
  merged <- FALSE # flag indicating components have been merged
  while (!merged) {
    fit <- mich_tree_cpp(y, edges, L, fit_intercept, fit_scale, tol, max_iter,
                         verbose = verbose & !L_auto,
                         log_pi_l, omega_l, mu_0, lambda_0,
                         b_bar_l, omega_bar_l, pi_bar_l, log_pi_bar_l)
    
    merged <- TRUE
    refit <- TRUE
    
    # identify components to merge
    if (L > 1) {
      # only merge columns with credible sets with length less than detect
      cred_sets <- apply(pi_bar_l, 2, cred_set, level=merge_level, simplify=FALSE)
      keep <- sapply(cred_sets, length) <= detect
      keep_mat <- matrix(FALSE, ncol = L, nrow = L)
      keep_mat[keep, keep] <- TRUE
      diag(keep_mat) <- FALSE

      # compute pairwise merge probabilities
      merge_prob_mat <- t(pi_bar_l) %*% pi_bar_l
      diag(merge_prob_mat) <- 0

      merge_residual <- fit$r_bar
      merge_delta <- fit$delta

      # merge mean components with pairwise merge probabilities > merge_prob
      while (L > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
        merged <- FALSE
        L <- L - 1

        # identify components with largest pairwise merge probabilities
        merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]

        # keep probabilities of component with largest posterior probability
        if (max(pi_bar_l[,merge_dex[2]]) > max(pi_bar_l[,merge_dex[1]])) {
          pi_bar_l[,merge_dex[1]] <- pi_bar_l[,merge_dex[2]]
          log_pi_bar_l[,merge_dex[1]] <- log_pi_bar_l[,merge_dex[2]]
          merge_prob_mat[merge_dex[1], ] <- merge_prob_mat[merge_dex[2], ]
          merge_prob_mat[, merge_dex[1]] <- merge_prob_mat[, merge_dex[2]]
        }

        # store merged parameters
        b_bar_l[,merge_dex[1]] <- b_bar_l[,merge_dex[1]] + b_bar_l[,merge_dex[2]]

        # drop merged component
        merge_prob_mat <- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
        keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
        pi_bar_l <- pi_bar_l[,-merge_dex[2], drop=FALSE]
        log_pi_bar_l <- log_pi_bar_l[,-merge_dex[2], drop=FALSE]
        b_bar_l <- b_bar_l[,-merge_dex[2], drop=FALSE]
        omega_bar_l <- omega_bar_l[,-merge_dex[2], drop=FALSE]
        if (L_auto) {
          log_pi_l <- log_pi_l[,-merge_dex[2], drop=FALSE]
        }
      }
    }
  }

  # if components were merged out use auto procedure to increase to desired LKJ
  merge_flag <- (L < L_max & !L_auto)
  merge_elbo <- -Inf
  if (merge_flag) increment <- 1 # ensures we don't overshoot number of components

  #### auto procedure with single component ####
  # 1. Increase L if ELBO increases, set restart == TRUE and increase Lagain
  # 2. If ELBO decreases, decrease counter and fit another component
  # 3. If ELBO decreases, counter == 0, and restart == TRUE, reset components
  #    to null model (except for this with concentrated probs) and refit
  # 4. If ELBO decreases, counter == 0, restart == FALSE, and merge_counter > 0,
  #    set fit to best model so far and merge components if no merges return
  #    model, otherwise decrease merge_count and go back to 1
  # 5. If merge_count = 0 return fit

  last_restart <- ifelse(restart, 2, Inf)

  if (L_auto | merge_flag) {

    refit <- L > 0
    counter <- n_search # number of searches after max elbo
    elbo <- fit$elbo[length(fit$elbo)] # current value of elbo
    elbo_new <- elbo

    if (verbose) print(paste0("L = ", L,": ELBO = ", elbo_new, "; Counter: ", counter))

    # continue search until n_search exhausted or max components exceeded
    while (L < L_max) {
      if (L_auto | L < L_max) {
        # increment dimension of parameters
        L <- L + increment
        if (L > 1 & L_auto) {
          log_pi_l <- cbind(matrix(log_pi_l[,1], nrow = n_node, ncol = increment), log_pi_l)
        }
        pi_bar_l <- cbind(matrix(1/n_node, nrow = n_node, ncol = increment), pi_bar_l)
        log_pi_bar_l <- cbind(matrix(0.0, nrow = n_node, ncol = increment), log_pi_bar_l)
        b_bar_l <- cbind(matrix(0.0, nrow = n_node, ncol = increment), b_bar_l)
        omega_bar_l <- cbind(matrix(1.0, nrow = n_node, ncol = increment), omega_bar_l)
      }

      # fit incremented model
      fit_new <- mich_tree_cpp(y, edges, L, fit_intercept, fit_scale, tol, max_iter,
                               verbose = FALSE, log_pi_l, omega_l, mu_0, lambda_0,
                               b_bar_l, omega_bar_l, pi_bar_l, log_pi_bar_l)

      # test if model improved or merge/restart
      if (!refit) refit <- TRUE
      elbo_new <- fit_new$elbo[length(fit_new$elbo)]
      if (verbose) print(paste0("L = ", L,": ELBO = ", elbo_new, "; Counter: ", counter))

      if (elbo_new > elbo | merge_flag) {
        elbo <- elbo_new
        fit <- fit_new
        counter <- n_search
        if (last_restart < Inf) restart <- TRUE
      } else if (counter == 0 & restart) {
        if (verbose) print(paste0("Restarting at L = ", L))
        restart <- FALSE
        counter <- n_search

        if (L_auto & L > last_restart) {
          last_restart <- L

          # reorder by longest blocks first
          chp <- apply(pi_bar_l, 2, which.max)
          chp_order <- order(chp)
          chp <- chp[chp_order]

          pi_bar_l <- pi_bar_l[,chp_order, drop = FALSE]
          log_pi_bar_l <- log_pi_bar_l[,chp_order, drop = FALSE]
          b_bar_l <- b_bar_l[,chp_order, drop = FALSE]
          omega_bar_l <- omega_bar_l[,chp_order, drop = FALSE]

          diff_order <- order(diff(c(chp, n_node)))

          pi_bar_l <- pi_bar_l[,diff_order, drop = FALSE]
          log_pi_bar_l <- log_pi_bar_l[,diff_order, drop = FALSE]
          b_bar_l <- b_bar_l[,diff_order, drop = FALSE]
          omega_bar_l <- omega_bar_l[,diff_order, drop = FALSE]

          # identify components with max prob > keep_level
          keep <- apply(pi_bar_l, 2, max) > keep_level

          # reset components with max prob < keep_level to null components
          if (sum(keep) < L) {
            keep_inc <- max(sum(!keep) - increment, 0)
            L <- sum(keep) + keep_inc

            log_pi_l <- sapply(1:L, function(i) log_pi_l[,1, drop = FALSE])
            pi_bar_l <- cbind(matrix(1/n_node, nrow = n_node, ncol = keep_inc), pi_bar_l[,keep, drop = FALSE])
            log_pi_bar_l <- cbind(matrix(0.0, nrow = n_node, ncol = keep_inc), log_pi_bar_l[,keep, drop = FALSE])
            b_bar_l <- cbind(matrix(0.0, nrow = n_node, ncol = keep_inc), b_bar_l[,keep, drop = FALSE])
            omega_bar_l <- cbind(matrix(1.0, nrow = n_node, ncol = keep_inc), omega_bar_l[,keep, drop = FALSE])
          }
        }

      } else if (counter == 0 | L == L_max) {
        # merging components
        counter <- n_search # reset counter in case searching continues
        merge_counter <- merge_counter - 1
        no_merges <- TRUE  # flag indicating no components were merged
        merged <- FALSE # flag indicating components have been merged
        if (verbose) print(paste0("Merging. Merge Counter: ", merge_counter))

        # set fit to best model so far
        L <- fit$L

        if (L >= 1) {
          log_pi_l <- log_pi_l[, 1:L, drop = FALSE]
          pi_bar_l <- fit$pi_bar_l
          log_pi_bar_l <- log(pi_bar_l)
          b_bar_l <- fit$b_bar_l
          omega_bar_l <- fit$omega_bar_l
        }
        
        while (!merged & L > 0) {
          fit <- mich_tree_cpp(y, edges, L, fit_intercept, fit_scale, tol, max_iter,
                               verbose = FALSE, log_pi_l, omega_l, mu_0, lambda_0,
                               b_bar_l, omega_bar_l, pi_bar_l, log_pi_bar_l)

          merged <- TRUE

          # identify components to merge
          if (L > 1) {
            # only merge columns with credible sets with length less than detect
            cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
            keep <- sapply(cred_sets, length) <= detect
            keep_mat <- matrix(FALSE, ncol = L, nrow = L)
            keep_mat[keep, keep] <- TRUE
            diag(keep_mat) <- FALSE

            # compute pairwise merge probabilities
            merge_prob_mat <- t(pi_bar_l) %*% pi_bar_l
            diag(merge_prob_mat) <- 0

            # merge mean components with pairwise merge probabilities > merge_prob
            while (L > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
              merged <- FALSE
              L <- L - 1

              # identify components with largest pairwise merge probabilities
              merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]

              # keep probabilities of component with largest posterior probability
              if (max(pi_bar_l[,merge_dex[2]]) > max(pi_bar_l[,merge_dex[1]])) {
                pi_bar_l[,merge_dex[1]] <- pi_bar_l[,merge_dex[2]]
                log_pi_bar_l[,merge_dex[1]] <- log_pi_bar_l[,merge_dex[2]]
                merge_prob_mat[merge_dex[1], ] <- merge_prob_mat[merge_dex[2], ]
                merge_prob_mat[, merge_dex[1]] <- merge_prob_mat[, merge_dex[2]]
              }

              # store merged parameters
              b_bar_l[,merge_dex[1]] <- b_bar_l[,merge_dex[1]] + b_bar_l[,merge_dex[2]]

              # drop merged component
              merge_prob_mat <- merge_prob_mat[-merge_dex[2], -merge_dex[2]]
              keep_mat <- keep_mat[-merge_dex[2], -merge_dex[2]]
              pi_bar_l <- pi_bar_l[,-merge_dex[2], drop=FALSE]
              log_pi_bar_l <- log_pi_bar_l[,-merge_dex[2], drop=FALSE]
              b_bar_l <- b_bar_l[,-merge_dex[2], drop=FALSE]
              omega_bar_l <- omega_bar_l[,-merge_dex[2], drop=FALSE]
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

        if (no_merges | merge_counter == 0) {
          fit <- fit_merge
          break
        }
      } else {
        counter <- counter - 1
      }
    }
  }

  #### return model ####
  class(fit) <- "mich.tree.fit"
  
  # rescale data to original units
  if (standardize) {
    fit$y <- fit$y * scale + center
    fit$r_bar <- fit$r_bar * scale
    fit$delta <- fit$delta * scale^2
    fit$mu_0 <- fit$mu_0 * scale + center
    fit$mu_bar <- fit$mu_bar * scale + center
    fit$lambda_0 <- fit$lambda_0 / scale^2

    if (fit$L > 0) {
      fit$b_bar_l <- fit$b_bar_l * scale
      fit$omega_bar_l <- fit$omega_bar_l / scale^2
      fit$mu_bar_l <- fit$mu_bar_l * scale
      fit$mu2_bar_l <- fit$mu2_bar_l * scale^2
    }
  }
  
  return(fit)
}

#' Single Change-Point Posterior Credible Set
#'
#' @description
#' The function `cred_set()` takes a length T vector of posterior change-point
#' location probabilities `prob` and a coverage level `level`, and returns the
#' smallest set of indices `s` such that `prob[s] > level`.
#'
#' @param prob  A numeric vector. A vector of posterior probabilities for the
#'   location of the change-point.
#' @param level A scalar. A single number in (0,1) that gives the lower bound
#'   for the probability that the credible set contains a change-point.
#'
#' @return A vector. A Level `level` posterior credible set for the location of
#' a single change-point.
#'
cred_set <- function(prob, level) {
  order(prob, decreasing = TRUE)[1:which.max(cumsum(prob[order(prob, decreasing = TRUE)]) > level)]
}

#' MICH Posterior Credible Sets
#'
#' @description
#' The function `mich_sets()` takes a T x N matrix of posterior change-point
#' location probabilities `probs`, a coverage level `level`, and a max set
#' length, and for column `i` of `probs` returns the MAP estimator
#' `which.max(probs[,i])` and the smallest set of indices `s_i` such that
#' `probs[s_i,i] > level` if `length(s_i) < max_length`.
#'
#' @param probs A numeric Matrix. A T x N matrix of posterior
#'   probabilities for the location of the change-points.
#' @param level A scalar. A single number in (0,1) that gives the lower bound
#'   for the probability that each credible set contains a change-point.
#' @param max_length A positive scalar. Detection threshold, if a credible set
#'   contains more that `max_length` indices, then no change is detected. Set
#'   equal to `log(T)^1.5` by default (see Section 2.5 of Berlind,
#'   Cappello, and Madrid Padilla (2025)).
#'
#' @return A list. MAP estimator of each change-point and corresponding credible
#' set.
#'
#' @export
#'
mich_sets <- function(probs, max_length = log(nrow(probs))^1.5, level = 0.9) {
  
  # initialize credible sets and changepoints
  cs <- list()
  est_cp <- numeric(0)
  
  cred_sets <- apply(probs, 2, cred_set, level = level, simplify = FALSE)
  
  # drop probs/sets that are longer than max_length
  keep <- which(sapply(cred_sets, length) <= max_length)
  if (length(keep) > 0) {
    est_cp <- apply(probs[, keep, drop = FALSE], 2, which.max)
    cs <- lapply(keep, function(i) cred_sets[[i]])
    # order change-points
    cs <- lapply(order(est_cp), function(i) cs[[i]])
    est_cp <- est_cp[order(est_cp)]
  }
  return(list(cp = est_cp, sets = lapply(cs, function(set) set[order(set)])))
}