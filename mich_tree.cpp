#include <Rcpp.h>
#include <unordered_map>
#include <vector>
using namespace Rcpp;

// Basic tree node structure
struct TreeNode {
  int id;
  double* r_bar;
  double* delta;
  double* scale;
  double mu, mu2;
  TreeNode* left = nullptr;
  TreeNode* right = nullptr;
  TreeNode* parent = nullptr;
  TreeNode* sibling = nullptr;

  TreeNode(int id_) : id(id_) {}
};

// Build the tree from an edge matrix
TreeNode* buildTree(NumericVector r_bar,
                    NumericVector delta,
                    double* lambda,
                    IntegerMatrix edges) {

  int n_leaf = r_bar.length();
  std::unordered_map<int, TreeNode*> nodes;
  std::unordered_map<int, std::vector<int>> children;

  for (int i = 0; i < edges.nrow(); i++) {
    int parent = edges(i,0) - n_leaf - 1;
    int child = edges(i,1) - n_leaf - 1;

    children[parent].push_back(child);

    if (!nodes.count(parent)) nodes[parent] = new TreeNode(parent);
    if (!nodes.count(child))  nodes[child]  = new TreeNode(child);
    if (nodes[parent]->left == nullptr) {
      nodes[parent]->left = nodes[child];
      nodes[child]->parent = nodes[parent];
    } else {
      nodes[parent]->right = nodes[child];
      nodes[child]->parent = nodes[parent];
      nodes[parent]->left->sibling = nodes[parent]->right;
      nodes[parent]->right->sibling = nodes[parent]->left;
    }
  }

  // Find root and assign leaf data
  int root;
  for (auto& kv : nodes) {
    TreeNode* node = kv.second;
    if (node->left == nullptr && node->right == nullptr) {
      node->id += n_leaf; // leaf ids go from 0:(n_leaf - 1)
      node->r_bar = &r_bar[node->id];
      node->delta = &delta[node->id];
      node->scale = lambda;
    } else if (node->parent == nullptr) {
      root = node->id;
    }
  }

  return nodes[root];
}

void tree_smcp(TreeNode* node,
               NumericVector b_bar,
               NumericVector omega_bar,
               NumericVector pi_bar,
               NumericVector log_pi_bar,
               NumericVector log_pi,
               double omega_0,
               int n_leaf,
               double* log_pi_max) {

  double b, omega, log_prob;

  if (node->left->left == nullptr && node->left->right == nullptr) {
    b_bar[node->id] += *(node->left->r_bar) * *(node->left->scale);
    omega_bar[node->id] += *(node->left->scale);
  } else {
    tree_smcp(node->left, b_bar, omega_bar, pi_bar, log_pi_bar, log_pi, omega_0, n_leaf, log_pi_max);

    b = b_bar[node->left->id];
    omega = omega_bar[node->left->id];

    b_bar[node->id] += b;
    omega_bar[node->id] += omega;

    omega += omega_0;
    omega_bar[node->left->id] = omega;
    b_bar[node->left->id] = b / omega;
    log_prob = log_pi[node->left->id] + 0.5 * (b * b / omega - std::log(omega));
    log_pi_bar[node->left->id] = log_prob;
    *log_pi_max = std::max(*log_pi_max, log_prob);
  }

  if (node->right->left == nullptr && node->right->right == nullptr) {
    b_bar[node->id] += *(node->right->r_bar) * *(node->right->scale);
    omega_bar[node->id] += *(node->right->scale);
  } else {
    tree_smcp(node->right, b_bar, omega_bar, pi_bar, log_pi_bar, log_pi, omega_0, n_leaf, log_pi_max);

    b = b_bar[node->right->id];
    omega = omega_bar[node->right->id];

    b_bar[node->id] += b;
    omega_bar[node->id] += omega;

    omega += omega_0;
    omega_bar[node->right->id] = omega;
    b_bar[node->right->id] = b / omega;
    log_prob = log_pi[node->right->id] + 0.5 * (b * b / omega - std::log(omega));
    log_pi_bar[node->right->id] = log_prob;
    *log_pi_max = std::max(*log_pi_max, log_prob);
  }

  if (node->parent == nullptr) {
    b = b_bar[node->id];
    omega = omega_0 + omega_bar[node->id];

    omega_bar[node->id] = omega;
    b_bar[node->id] = b / omega;
    log_prob = log_pi[node->id] + 0.5 * (b * b / omega - std::log(omega));
    log_pi_bar[node->id] = log_prob;
    *log_pi_max = std::max(*log_pi_max, log_prob);

    double tot = 0;
    for (int i = 0; i < pi_bar.length(); i++) {
      log_pi_bar[i] -= *log_pi_max;
      if (log_pi[i] > R_NegInf) pi_bar[i] = std::exp(log_pi_bar[i]);
      else pi_bar[i] = 0.0;
      tot += pi_bar[i];
    }
    pi_bar = pi_bar / tot;
  }
  return;
}

void mu_bar_update(TreeNode* node,
                   NumericMatrix mu_bar,
                   NumericMatrix mu2_bar,
                   int l,
                   NumericVector b_bar,
                   NumericVector omega_bar,
                   NumericVector pi_bar) {

  if (node->left == nullptr && node->right == nullptr) {
    mu_bar(node->id, l) = node->mu;
    mu2_bar(node->id, l) = node->mu2;
    //*(node->r_bar) -= node->mu;
  } else {
    node->mu += b_bar[node->id] * pi_bar[node->id];
    node->mu2 += (1 / omega_bar[node->id] + b_bar[node->id] * b_bar[node->id]) * pi_bar[node->id];
    node->left->mu = node->mu;
    node->right->mu = node->mu;
    node->left->mu2 = node->mu2;
    node->right->mu2 = node->mu2;
    mu_bar_update(node->left, mu_bar, mu2_bar, l, b_bar, omega_bar, pi_bar);
    mu_bar_update(node->right, mu_bar, mu2_bar, l, b_bar, omega_bar, pi_bar);
  }
  return;
}


double elbo_fn(double lambda_0,
               NumericVector r_bar,
               NumericVector delta,
               NumericMatrix b_bar_l,
               NumericMatrix omega_bar_l,
               NumericMatrix pi_bar_l,
               NumericMatrix log_pi_bar_l,
               NumericMatrix log_pi_l,
               double omega_l,
               double log_omega_l)  {
  int L = pi_bar_l.ncol();
  int n_leaf = r_bar.length();
  int n_node = n_leaf - 1;

  // calculate E[log p(y)]
  double elbo = 0.5 * n_leaf * std::log(lambda_0) - 0.5 * lambda_0 * (r_bar[n_node] * r_bar[n_node] + delta[n_node]);
  for (int i = 0; i < n_node; i++) {
    elbo += -0.5 * lambda_0 * (r_bar[i] * r_bar[i] + delta[i]);
    for (int l = 0; l < L; l++) {
      if (pi_bar_l(i, l) > 1e-20) elbo += pi_bar_l(i, l) * (log_pi_l(i, l) - log_pi_bar_l(i, l));
      elbo += 0.5 * pi_bar_l(i, l) * (log_omega_l - std::log(omega_bar_l(i, l)));
      elbo += 0.5 * pi_bar_l(i, l) * (1 - omega_l * (b_bar_l(i, l) * b_bar_l(i, l) + 1 / omega_bar_l(i, l)));
    }
  }
  return elbo;
}

NumericVector tree_vb(TreeNode* tree,
                      NumericVector r_bar,
                      NumericVector delta,
                      bool fit_intercept,
                      bool fit_scale,
                      double tol,
                      double max_iter,
                      bool verbose,
                      double* mu_0,
                      double* lambda_0,
                      double omega_l,
                      NumericMatrix log_pi_l,
                      NumericMatrix mu_bar_l,
                      NumericMatrix mu2_bar_l,
                      NumericMatrix b_bar_l,
                      NumericMatrix omega_bar_l,
                      NumericMatrix pi_bar_l,
                      NumericMatrix log_pi_bar_l) {

  int n_leaf = r_bar.length(), iter = 0, L = log_pi_bar_l.ncol();
  double log_pi_max, log_omega_l = std::log(omega_l);
  
  // initialize ELBO
  std::vector<double> elbo;
  elbo.reserve(std::min(max_iter, 1e6));
  elbo.push_back(R_NegInf);

  NumericVector b_bar(n_leaf - 1), omega_bar(n_leaf - 1), pi_bar(n_leaf - 1), log_pi_bar(n_leaf - 1);

  while (iter < max_iter) {
    for (int l=0; l < L; l++) {
      r_bar += mu_bar_l(_, l);
      delta += (mu_bar_l(_, l) * mu_bar_l(_, l) - mu2_bar_l(_, l));

      b_bar.fill(0.0); omega_bar.fill(0.0); pi_bar.fill(0.0); log_pi_bar.fill(0.0);
      log_pi_max = R_NegInf;

      tree_smcp(tree, b_bar, omega_bar, pi_bar, log_pi_bar, log_pi_l(_, l), omega_l, n_leaf, &log_pi_max);

      // merge component into mu_0 if root is active
      if (fit_intercept) {
        bool merge_root = true;
        for (int i=1; i < n_leaf - 1; i++) {
          if (pi_bar[i] > pi_bar[0]) {
            merge_root = false;
            break;
          }
        }
        if (merge_root) {
          *mu_0 += pi_bar[0] * b_bar[0];
          r_bar += Rcpp::rep(-*mu_0, n_leaf);
          b_bar.fill(0.0); 
          omega_bar.fill(1.0); 
          pi_bar.fill(1.0 / (n_leaf - 1)); 
          log_pi_bar.fill(0.0);
        }
      }
      
      b_bar_l(_, l) = b_bar;
      omega_bar_l(_, l) = omega_bar;
      pi_bar_l(_, l) = pi_bar;
      log_pi_bar_l(_, l) = log_pi_bar;

      tree->mu = 0.0;
      tree->mu2 = 0.0;
      mu_bar_update(tree, mu_bar_l, mu2_bar_l, l, b_bar, omega_bar, pi_bar);

      r_bar += -mu_bar_l(_, l);
      delta += (mu2_bar_l(_, l) - mu_bar_l(_, l) * mu_bar_l(_, l));
    }

    if (fit_intercept) {
      r_bar += Rcpp::rep(*mu_0, n_leaf);
      *mu_0 = Rcpp::mean(r_bar);
      r_bar += Rcpp::rep(-*mu_0, n_leaf);
    }
    
    if (fit_scale) {
      *lambda_0 = n_leaf / Rcpp::sum((Rcpp::pow(r_bar, 2) + delta));
    }

    iter++;
    elbo.push_back(0.0);
    Rcout << "lambda_0:" << *lambda_0 << ";\n";
    Rcout << "r_bar:" << r_bar[0] << ";\n";
    Rcout << "b_bar_l:" << b_bar_l(0,0) << ";\n";
    Rcout << "omega_bar_l:" << omega_bar_l(0,0) << ";\n";
    Rcout << "pi_bar_l:" << pi_bar_l(0,0) << ";\n";
    Rcout << "pi_bar_l:" << pi_bar_l(0,0) << ";\n";
    Rcout << "log_pi_bar_l:" << log_pi_bar_l(0,0) << ";\n";
    
    elbo[iter] = elbo_fn(*lambda_0, r_bar, delta, b_bar_l, omega_bar_l, pi_bar_l,
                         log_pi_bar_l, log_pi_l, omega_l, log_omega_l);
    
    if (std::isnan(elbo[iter])) throw std::runtime_error("NaN in elbo");
    if (verbose && (iter % 5000 == 0)) Rcout << "Iteration: " << iter << "; elbo: " << elbo[iter] << ";\n";
    // terminate if relative increase in elbo falls below tol
    if (std::abs((elbo[iter] - elbo[iter - 1]) / elbo[iter - 1]) < tol) break;
  }
  return wrap(elbo);
}

void sibling_merge(TreeNode* node,
                   IntegerVector active,
                   NumericMatrix b_bar,
                   NumericMatrix omega_bar,
                   NumericMatrix pi_bar,
                   NumericMatrix log_pi_bar) {

  if (node->left == nullptr && node->right == nullptr) return;

  sibling_merge(node->left, active, b_bar, omega_bar, pi_bar, log_pi_bar);
  sibling_merge(node->right, active, b_bar, omega_bar, pi_bar, log_pi_bar);

  bool left = false, right = false;
  int left_dex, right_dex;

  for (int i = 0; i < active.length(); i++) {
    if (node->left->id == active[i]) {
      left = true;
      left_dex = i;
    } else if (node->right->id == active[i]) {
      right = true;
      right_dex = i;
    }
  }

  if (left && right) {
    double diff = std::abs(b_bar(node->left->id, left_dex) - b_bar(node->right->id, right_dex));
    double diff_se = std::sqrt(1 / omega_bar(node->left->id, left_dex) + 1 / omega_bar(node->right->id, right_dex));
    if (diff - 1.96 * diff_se < 0 && diff + 1.96 * diff_se > 0) {
      // merge children both into parent
      active[right_dex] = node->id;
      active[left_dex] = -1;

      double b_tmp = 0.5 * (b_bar(node->right->id, right_dex) + b_bar(node->left->id, left_dex));

      if (pi_bar(node->right->id, right_dex) > pi_bar(node->left->id, left_dex)) {
        b_bar(node->right->id, right_dex) = b_bar(node->id, right_dex);
        b_bar(node->id, right_dex) = b_tmp;

        double pi_tmp = pi_bar(node->right->id, right_dex);
        pi_bar(node->right->id, right_dex) = pi_bar(node->id, right_dex);
        pi_bar(node->id, right_dex) = pi_tmp;

        double log_pi_tmp = log_pi_bar(node->right->id, right_dex);
        log_pi_bar(node->right->id, right_dex) = log_pi_bar(node->id, right_dex);
        log_pi_bar(node->id, right_dex) = log_pi_tmp;

        for (int i=0; i < b_bar.nrow(); i++) {
          b_bar(i, left_dex) = 0.0;
          pi_bar(i, left_dex) = 1.0 / b_bar.nrow();
          log_pi_bar(i, left_dex) = 0.0;
        }
      } else {
        b_bar(node->left->id, left_dex) = b_bar(node->id, left_dex);
        b_bar(node->id, left_dex) = b_tmp;

        double pi_tmp = pi_bar(node->left->id, left_dex);
        pi_bar(node->left->id, left_dex) = pi_bar(node->id, left_dex);
        pi_bar(node->id, left_dex) = pi_tmp;

        double log_pi_tmp = log_pi_bar(node->left->id, left_dex);
        log_pi_bar(node->left->id, left_dex) = log_pi_bar(node->id, left_dex);
        log_pi_bar(node->id, left_dex) = log_pi_tmp;

        for (int i=0; i < b_bar.nrow(); i++) {
          b_bar(i, right_dex) = 0.0;
          pi_bar(i, right_dex) = 1.0 / b_bar.nrow();
          log_pi_bar(i, right_dex) = 0.0;
        }
      }
    } else {
      for (int i=0; i < active.length(); i++) {
        if (node->sibling->id == active[i]) {
          if (std::abs(b_bar(node->right->id, right_dex) - b_bar(node->sibling->id, i)) < std::abs(b_bar(node->left->id, left_dex) - b_bar(node->sibling->id, i))) {
            // only merge right child into parent
            active[right_dex] = node->id;

            double b_tmp = b_bar(node->right->id, right_dex);
            b_bar(node->right->id, right_dex) = b_bar(node->id, right_dex);
            b_bar(node->id, right_dex) = b_tmp;
            b_bar(node->left->id, left_dex) -= b_tmp;

            double pi_tmp = pi_bar(node->right->id, right_dex);
            pi_bar(node->right->id, right_dex) = pi_bar(node->id, right_dex);
            pi_bar(node->id, right_dex) = pi_tmp;

            double log_pi_tmp = log_pi_bar(node->right->id, right_dex);
            log_pi_bar(node->right->id, right_dex) = log_pi_bar(node->id, right_dex);
            log_pi_bar(node->id, right_dex) = log_pi_tmp;
            return;
          }
        }
      }

      // only merge left child into parent
      active[left_dex] = node->id;

      double b_tmp = b_bar(node->left->id, left_dex);
      b_bar(node->left->id, left_dex) = b_bar(node->id, left_dex);
      b_bar(node->id, left_dex) = b_tmp;
      b_bar(node->right->id, right_dex) -= b_tmp;

      double pi_tmp = pi_bar(node->left->id, left_dex);
      pi_bar(node->left->id, left_dex) = pi_bar(node->id, left_dex);
      pi_bar(node->id, left_dex) = pi_tmp;

      double log_pi_tmp = log_pi_bar(node->left->id, left_dex);
      log_pi_bar(node->left->id, left_dex) = log_pi_bar(node->id, left_dex);
      log_pi_bar(node->id, left_dex) = log_pi_tmp;
    }
  }

  return;
}

void merge_reset(TreeNode* node,
                   IntegerVector active,
                   LogicalVector merged_up,
                   NumericMatrix b_bar,
                   NumericMatrix omega_bar,
                   NumericMatrix pi_bar,
                   NumericMatrix log_pi_bar) {

  if (node->left == nullptr && node->right == nullptr) return;

  merge_reset(node->left, active, merged_up, b_bar, omega_bar, pi_bar, log_pi_bar);
  merge_reset(node->right, active, merged_up, b_bar, omega_bar, pi_bar, log_pi_bar);

  bool left = false, right = false;
  int left_dex, right_dex;

  for (int i = 0; i < active.length(); i++) {
    if (node->left->id == active[i]) {
      left = true;
      left_dex = i;
    } else if (node->right->id == active[i]) {
      right = true;
      right_dex = i;
    }
  }

  if (left && right) {
    if (merged_up[left_dex]) {
      // merge right child into parent and reset left clade
      active[right_dex] = node->id;
      merged_up[right_dex] = true;

      // for (int i=0; i < b_bar.nrow(); i++) {
      //   if (i != node->id) {
      //     b_bar(node->id, right_dex) += b_bar(i, right_dex) * pi_bar(i, right_dex) / pi_bar(node->right->id, right_dex);
      //     b_bar(i, right_dex) = 0.0;
      //   }
      // }

      double b_tmp = b_bar(node->right->id, right_dex);
      b_bar(node->right->id, right_dex) = b_bar(node->id, right_dex);
      b_bar(node->id, right_dex) = b_tmp;
      
      double pi_tmp = pi_bar(node->right->id, right_dex);
      pi_bar(node->right->id, right_dex) = pi_bar(node->id, right_dex);
      pi_bar(node->id, right_dex) = pi_tmp;

      double log_pi_tmp = log_pi_bar(node->right->id, right_dex);
      log_pi_bar(node->right->id, right_dex) = log_pi_bar(node->id, right_dex);
      log_pi_bar(node->id, right_dex) = log_pi_tmp;

      for (int i = 0; i < active.length(); i++) {
        if (node->left->id == active[i] || node->left->left->id == active[i] || node->left->right->id == active[i]) {
          active[i] = -1;
          b_bar(_, i) = NumericVector(b_bar.nrow(), 0.0);
          pi_bar(_, i) = NumericVector(pi_bar.nrow(), 1.0 / pi_bar.nrow());
          log_pi_bar(_, i) =  NumericVector(log_pi_bar.nrow(), 0.0);
        }
      }
    } else if (merged_up[right_dex]) {
      // merge left child into parent and reset right clade
      active[left_dex] = node->id;
      merged_up[left_dex] = true;

      // for (int i=0; i < b_bar.nrow(); i++) {
      //   if (i != node->id) {
      //     b_bar(node->id, left_dex) += b_bar(i, left_dex) * pi_bar(i, left_dex) / pi_bar(node->right->id, left_dex);
      //     b_bar(i, left_dex) = 0.0;
      //   }
      // }
      
      double b_tmp = b_bar(node->left->id, left_dex);
      b_bar(node->left->id, left_dex) = b_bar(node->id, left_dex);
      b_bar(node->id, left_dex) = b_tmp;
      
      double pi_tmp = pi_bar(node->left->id, left_dex);
      pi_bar(node->left->id, left_dex) = pi_bar(node->id, left_dex);
      pi_bar(node->id, left_dex) = pi_tmp;

      double log_pi_tmp = log_pi_bar(node->left->id, left_dex);
      log_pi_bar(node->left->id, left_dex) = log_pi_bar(node->id, left_dex);
      log_pi_bar(node->id, left_dex) = log_pi_tmp;

      for (int i = 0; i < active.length(); i++) {
        if (node->right->id == active[i] || node->right->left->id == active[i] || node->right->right->id == active[i]) {
          active[i] = -1;
          b_bar(_, i) = NumericVector(b_bar.nrow(), 0.0);
          pi_bar(_, i) = NumericVector(pi_bar.nrow(), 1.0 / pi_bar.nrow());
          log_pi_bar(_, i) = NumericVector(log_pi_bar.nrow(), 0.0);
        }
      }
    } else {
      // merge right child into parent
      active[right_dex] = node->id;
      merged_up[right_dex] = true;

      double b_tmp = b_bar(node->right->id, right_dex);
      b_bar(node->right->id, right_dex) = b_bar(node->id, right_dex);
      b_bar(node->id, right_dex) = b_tmp;
      b_bar(node->left->id, left_dex) -= b_tmp;

      double pi_tmp = pi_bar(node->right->id, right_dex);
      pi_bar(node->right->id, right_dex) = pi_bar(node->id, right_dex);
      pi_bar(node->id, right_dex) = pi_tmp;

      double log_pi_tmp = log_pi_bar(node->right->id, right_dex);
      log_pi_bar(node->right->id, right_dex) = log_pi_bar(node->id, right_dex);
      log_pi_bar(node->id, right_dex) = log_pi_tmp;
    }
  }
  return;
}

// [[Rcpp::export]]
List mich_tree_cpp(NumericVector y,
               IntegerMatrix edges,
               int L,
               bool fit_intercept,
               bool fit_scale,
               double tol,
               double max_iter,
               bool verbose,
               NumericMatrix log_pi_l,
               double omega_l,
               NumericVector mu_0_vec,
               NumericVector lambda_0_vec,
               NumericMatrix b_bar_l,
               NumericMatrix omega_bar_l,
               NumericMatrix pi_bar_l,
               NumericMatrix log_pi_bar_l) {

  int n_leaf = y.length();
  int n_node = n_leaf - 1;

  NumericVector r_bar = clone(y);
  NumericVector delta(n_leaf, 0.0);
  NumericMatrix mu_bar_l(n_leaf, L), mu2_bar_l(n_leaf, L);
  
  double mu_0 = mu_0_vec[0];
  double lambda_0 = lambda_0_vec[0];
  
  TreeNode* tree = buildTree(r_bar, delta, &lambda_0, edges);

  for (int l = 0; l < L; l++) {
    mu_bar_update(tree, mu_bar_l, mu2_bar_l, l, b_bar_l(_,l), omega_bar_l(_,l), pi_bar_l(_,l));
  }

  for (int i=0; i < n_leaf; i++) {
    r_bar[i] -= mu_0;
    for (int l=0; l < L; l++) {
      r_bar[i] -= mu_bar_l(i, l);
      delta[i] += mu2_bar_l(i, l) - mu_bar_l(i, l) * mu_bar_l(i, l);
    }
  }
  
  NumericVector elbo = tree_vb(tree, r_bar, delta, fit_intercept, fit_scale, tol,
                               max_iter, verbose, &mu_0, &lambda_0, omega_l,
                               log_pi_l, mu_bar_l, mu2_bar_l, b_bar_l,
                               omega_bar_l, pi_bar_l, log_pi_bar_l);
  
  IntegerVector active(L);
  LogicalVector merged_up(L, false);

  for (int l=0; l < L; l++) {
    double pi_max = 0.0;
    for (int i=0; i < n_node; i++) {
      if (pi_bar_l(i, l) > pi_max) {
        pi_max = pi_bar_l(i, l);
        active[l] = i;
      }
    }
  }
  
  merge_reset(tree, active, merged_up, b_bar_l, omega_bar_l, pi_bar_l, log_pi_bar_l);
  //sibling_merge(tree, active, b_bar_l, omega_bar_l, pi_bar_l, log_pi_bar_l);
  for (int l = 0; l < L; l++) {
    r_bar += mu_bar_l(_,l);
    delta += mu_bar_l(_,l) * mu_bar_l(_,l) - mu2_bar_l(_,l);
    mu_bar_update(tree, mu_bar_l, mu2_bar_l, l, b_bar_l(_,l), omega_bar_l(_,l), pi_bar_l(_,l));
    r_bar += -mu_bar_l(_,l);
    delta += mu2_bar_l(_,l) - mu_bar_l(_,l) * mu_bar_l(_,l);
  }

  NumericVector elbo_merged = tree_vb(tree, r_bar, delta, fit_intercept, fit_scale, tol,
                                      max_iter, verbose, &mu_0, &lambda_0, omega_l,
                                      log_pi_l, mu_bar_l, mu2_bar_l, b_bar_l,
                                      omega_bar_l, pi_bar_l, log_pi_bar_l);

  for (int i = 0; i < elbo_merged.length(); i++) {
    elbo.push_back(elbo_merged[i]);
  }

  merged_up.fill(false);
  for (int l=0; l < L; l++) {
    double pi_max = 0.0;
    for (int i=0; i < n_node; i++) {
      if (pi_bar_l(i, l) > pi_max) {
        pi_max = pi_bar_l(i, l);
        active[l] = i;
      }
    }
  }

  merge_reset(tree, active, merged_up, b_bar_l, omega_bar_l, pi_bar_l, log_pi_bar_l);
  //sibling_merge(tree, active, b_bar_l, omega_bar_l, pi_bar_l, log_pi_bar_l);
  for (int l = 0; l < L; l++) {
    r_bar += mu_bar_l(_,l);
    delta += mu_bar_l(_,l) * mu_bar_l(_,l) - mu2_bar_l(_,l);
    mu_bar_update(tree, mu_bar_l, mu2_bar_l, l, b_bar_l(_,l), omega_bar_l(_,l), pi_bar_l(_,l));
    r_bar += -mu_bar_l(_,l);
    delta += mu2_bar_l(_,l) - mu_bar_l(_,l) * mu_bar_l(_,l);
  }

  mu_0_vec[0] = mu_0;
  lambda_0_vec[0] = lambda_0;
  
  return List::create(_["y"] = y,
                      _["r_bar"] = r_bar,
                      _["delta"] = delta,
                      _["mu_0"] = mu_0,
                      _["lambda_0"] = lambda_0,
                      _["mu_bar"] = Rcpp::rowSums(mu_bar_l) + Rcpp::rep(mu_0, n_leaf),
                      _["mu_bar_l"] = mu_bar_l,
                      _["mu2_bar_l"] = mu2_bar_l,
                      _["b_bar_l"] = b_bar_l,
                      _["omega_bar_l"] = omega_bar_l,
                      _["pi_bar_l"] = pi_bar_l,
                      _["log_pi_bar_l"] = log_pi_bar_l,
                      _["elbo"] = elbo,
                      _["L"] = L);
}

bool contains_int(IntegerVector vec, int val) {
  for (int i = 0; i < vec.length(); i++) {
    if (vec[i] == val) return true;
  }
  return false;
}