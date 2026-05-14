#include <Rcpp.h>
#include <unordered_map>
#include <vector>
using namespace Rcpp;

// Basic tree node structure
struct TreeNode {
  int id;
  double n_leaf;
  NumericVector QTy;
  NumericVector eigen_val;
  NumericVector mu;
  double muTmu;
  
  // graph navigation pointers
  TreeNode* left = nullptr;
  TreeNode* right = nullptr;
  TreeNode* parent = nullptr;
  TreeNode* sibling = nullptr;
  
  TreeNode(int id_, int size) : 
    id(id_),
    QTy(size),
    eigen_val(size),
    mu(size)
  {}
};

void leaf_count(TreeNode* node) {
  if (node->left == nullptr && node->right == nullptr) {
    node->n_leaf = 1.0; 
  } else {
    leaf_count(node->left);
    leaf_count(node->right);
    node->n_leaf = node->left->n_leaf + node->right->n_leaf;
  }
  return;
};

void eigen_update(
    TreeNode* node, 
    NumericVector lambda, 
    NumericVector Lambda_bar_log_det,
    NumericVector Lambda_bar_trace,
    NumericMatrix mean_wts,
    NumericMatrix sandwich_wts,
    double omega_0
  ) {
  
  int d =  lambda.length();
  NumericVector eigen_val;
  
  if (node->left == nullptr && node->right == nullptr) {
    eigen_val = omega_0 + lambda;
  } else {
    eigen_update(node->left, lambda, Lambda_bar_log_det, Lambda_bar_trace, mean_wts, sandwich_wts, omega_0);
    eigen_update(node->right, lambda, Lambda_bar_log_det, Lambda_bar_trace, mean_wts, sandwich_wts, omega_0);
    eigen_val = node->left->eigen_val + node->right->eigen_val - omega_0; // omega_0 gets double counted
  }
  
  node->eigen_val = eigen_val;
  
  Lambda_bar_log_det[node->id] = 0.0; 
  Lambda_bar_trace[node->id] = 0.0; 
  for (int i = 0; i < d; i++) {
    Lambda_bar_log_det[node->id] += std::log(eigen_val[i]);
    Lambda_bar_trace[node->id] += eigen_val[i];
  }
  
  mean_wts(node->id, _) = lambda / eigen_val; 
  sandwich_wts(node->id, _) = lambda * mean_wts(node->id, _); 
  return;
};

// Build the tree from an edge matrix
TreeNode* buildTree(
    IntegerMatrix edges,
    NumericVector lambda,
    double omega_0
  ) {
  
  int d = lambda.size();
  
  // create list of nodes and children
  std::unordered_map<int, TreeNode*> nodes;
  std::unordered_map<int, std::vector<int>> children;
  
  for (int i = 0; i < edges.nrow(); i++) {
    // extract parent and child node index
    int parent = edges(i,0) - 1;
    int child = edges(i,1) - 1;
    
    children[parent].push_back(child);
    
    // initialize parent if necessary
    if (!nodes.count(parent)) nodes[parent] = new TreeNode(parent, d);
    // initialize child if necessary
    if (!nodes.count(child))  nodes[child]  = new TreeNode(child, d);
    
    // build graph left to right
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
  
  // Find root
  int root;
  for (auto& kv : nodes) {
    TreeNode* node = kv.second;
    if (node->parent == nullptr) {
      root = node->id;
      break;
    }
  }
  
  // count leaf offspring for each node
  leaf_count(nodes[root]);
  
  return nodes[root];
}

void mu_bar_rec(
    TreeNode* node,
    NumericVector lambda,
    NumericMatrix mean_wts,
    NumericMatrix mu_bar,
    NumericVector muTmu_bar,
    NumericMatrix b_bar,
    NumericVector pi_bar
  ) {
  // add internal node contribution to mean of each leaf
  node->mu += b_bar(node->id, _) * pi_bar[node->id];
  for (int i = 0; i < mu_bar.ncol(); i++) {
    node->muTmu += (mean_wts(node->id, i) + lambda[i] * b_bar(node->id, i) * b_bar(node->id, i)) * pi_bar[node->id];
  }
  
  if (node->left == nullptr && node->right == nullptr) {
    mu_bar(node->id, _) = node->mu; // store mean of leaf
    muTmu_bar[node->id] = node->muTmu;
  } else {
    // init offspring mean and recursively update to leafs
    node->left->mu = clone(node->mu);
    node->left->muTmu = node->muTmu;
    mu_bar_rec(node->left, lambda, mean_wts, mu_bar, muTmu_bar, b_bar, pi_bar); 
    node->right->mu = clone(node->mu);
    node->right->muTmu = node->muTmu;
    mu_bar_rec(node->right, lambda, mean_wts, mu_bar, muTmu_bar, b_bar, pi_bar); 
  }
  return;
}

void tree_multi_smcp_rec(
    TreeNode* node,
    NumericMatrix QTy,
    NumericVector log_det,
    NumericMatrix mean_wts,
    NumericMatrix sandwich_wts,
    NumericMatrix QTb_bar,
    NumericVector pi_bar,
    NumericVector log_pi_bar,
    NumericVector log_pi,
    double* log_pi_max
  ) {
  double QTy_i;
  int d = QTy.ncol();
  
  if (node->left == nullptr && node->right == nullptr) {
    node->QTy = QTy(node->id, _);
  } else {
    tree_multi_smcp_rec(node->left, QTy, log_det, mean_wts, sandwich_wts, QTb_bar, pi_bar, log_pi_bar, log_pi, log_pi_max);
    tree_multi_smcp_rec(node->right, QTy, log_det, mean_wts, sandwich_wts, QTb_bar, pi_bar, log_pi_bar, log_pi, log_pi_max);
    node->QTy = node->left->QTy + node->right->QTy;
  }
  
  log_pi_bar[node->id] = log_pi[node->id] - 0.5 * log_det[node->id];
  
  for (int i = 0; i < d; i++) {
    QTy_i = node->QTy[i];
    log_pi_bar[node->id] += 0.5 * sandwich_wts(node->id, i) * QTy_i * QTy_i;
    QTb_bar(node->id, i) = mean_wts(node->id, i) * QTy_i;
  }
  
  *log_pi_max = std::max(*log_pi_max, log_pi_bar[node->id]);
  return;
}

List tree_multi_smcp(
    TreeNode* root,
    NumericVector lambda,
    NumericMatrix QTy,
    NumericVector log_det,
    NumericMatrix mean_wts,
    NumericMatrix sandwich_wts,
    NumericVector log_pi
  ) { 
  int n_leaf = QTy.nrow();
  int n_node = log_pi.length();
  int d = QTy.ncol();
  double log_pi_max = R_NegInf;
  NumericMatrix QTb_bar (n_node, d);
  NumericVector pi_bar (n_node, 0.0);
  NumericVector log_pi_bar (n_node, 0.0);
  NumericMatrix QTmu_bar (n_leaf, d);
  NumericVector muTLmu_bar (n_leaf, 0.0);
  
  tree_multi_smcp_rec(root, QTy, log_det, mean_wts, sandwich_wts, QTb_bar, pi_bar, log_pi_bar, log_pi, &log_pi_max);
  
  // rescale probabilities 
  double tot = 0;
  for (int i = 0; i < n_node; i++) {
    log_pi_bar[i] -= log_pi_max;
    if (log_pi[i] > R_NegInf) pi_bar[i] = std::exp(log_pi_bar[i]);
    else pi_bar[i] = 0.0;
    tot += pi_bar[i];
  }
  pi_bar = pi_bar / tot;
  
  // calculate mean params
  root->mu.fill(0.0); // zero out mean before calling recursive update
  root->muTmu = 0.0;
  mu_bar_rec(root, lambda, mean_wts, QTmu_bar, muTLmu_bar, QTb_bar, pi_bar);
  
  return List::create(
    _["QTb_bar"] = QTb_bar,
    _["QTmu_bar"] = QTmu_bar,
    _["muTLmu_bar"] = muTLmu_bar,
    _["pi_bar"] = pi_bar,
    _["log_pi_bar"] = log_pi_bar
  );
}

double tree_multi_elbo(
  NumericVector lambda,
  NumericMatrix QTr_bar,
  List post_params,
  NumericVector Lambda_bar_log_det,
  NumericVector Lambda_bar_trace,
  NumericMatrix log_pi_l,
  double omega,
  double d_log_omega
  ) {

  int n_leaf = QTr_bar.nrow();
  int n_node = Lambda_bar_log_det.length();
  int d = QTr_bar.ncol();
  int L = post_params.length();

  List post_params_l;
  NumericMatrix QTmu_bar_l (n_leaf, d);
  NumericVector muTLmu_bar_l (n_leaf);
  NumericVector pi_bar_l (n_node);
  NumericVector log_pi_bar_l (n_node);
  NumericMatrix QTb_bar_l (n_node, d);

  // E[log p]
  double elbo = 0.0;
  for (int j = 0; j < d; j++) {
    elbo += n_leaf * std::log(lambda[j]);
    for (int i = 0; i < n_leaf; i++) {
      elbo -= lambda[j] * QTr_bar(i, j) * QTr_bar(i, j);
    }
  }

  for (int l = 0; l < L; l++) {
    post_params_l = post_params[l];
    QTmu_bar_l = as<NumericMatrix>(post_params_l["QTmu_bar"]);
    muTLmu_bar_l = as<NumericVector>(post_params_l["muTLmu_bar"]);
    for (int i = 0; i < n_leaf; i++) {
      elbo -= muTLmu_bar_l[i];
      for (int j = 0; j < d; j++) {
        elbo += lambda[j] * QTmu_bar_l(i, j) * QTmu_bar_l(i, j);
      }
    }

    // KL(q || p)
    pi_bar_l = as<NumericVector>(post_params_l["pi_bar"]);
    log_pi_bar_l = as<NumericVector>(post_params_l["log_pi_bar"]);
    QTb_bar_l = as<NumericMatrix>(post_params_l["QTb_bar"]);
    for (int i = 0; i < n_node; i++) {
      if (pi_bar_l[i] > 1e-20) elbo += 2.0 * pi_bar_l[i] * (log_pi_l(i,l) - log_pi_bar_l[i]);
      elbo += pi_bar_l[i] * (Lambda_bar_log_det[i] - d_log_omega);
      elbo -= pi_bar_l[i] * (Lambda_bar_trace[i] / omega); // TODO cache divide by omega
      for (int j = 0; j < d; j++) {
        elbo -= omega * QTb_bar_l(i, j) * QTb_bar_l(i, j);
      }
    }
  }
  return elbo;
}

void merge_reset_rec(
    TreeNode* node,
    int n_node,
    int d,
    IntegerVector active,
    LogicalVector merged_up,
    List post_params
  ) {
  
  if (node->left == nullptr && node->right == nullptr) return;
  
  merge_reset_rec(node->left, n_node, d, active, merged_up, post_params);
  merge_reset_rec(node->right, n_node, d, active, merged_up, post_params);
  
  bool left = false, right = false;
  int left_dex, right_dex;
  int L = active.length();
  
  for (int l = 0; l < L; l++) {
    if (node->left->id == active[l]) {
      left = true;
      left_dex = l;
    } else if (node->right->id == active[l]) {
      right = true;
      right_dex = l;
    }
  }
  
  if (left && right) {
    if (merged_up[right_dex]) {
      // merge left child into parent and reset right clade
      active[left_dex] = node->id;
      merged_up[left_dex] = true;
      
      List post_params_l = post_params[left_dex];
      
      NumericMatrix QTb_bar_l = as<NumericMatrix>(post_params_l["QTb_bar"]);
      NumericVector QTb_tmp = QTb_bar_l(node->left->id, _);
      QTb_bar_l(node->left->id, _) = QTb_bar_l(node->id, _);
      QTb_bar_l(node->id, _) = QTb_tmp;
      post_params_l["QTb_bar"] = QTb_bar_l;
      
      NumericVector pi_bar_l = as<NumericVector>(post_params_l["pi_bar"]);
      double pi_tmp = pi_bar_l[node->left->id];
      pi_bar_l[node->left->id] = pi_bar_l[node->id];
      pi_bar_l[node->id] = pi_tmp;
      post_params_l["pi_bar"] = pi_bar_l;
      
      NumericVector log_pi_bar_l = as<NumericVector>(post_params_l["log_pi_bar"]);
      double log_pi_tmp = log_pi_bar_l[node->left->id];
      log_pi_bar_l[node->left->id] = log_pi_bar_l[node->id];
      log_pi_bar_l[node->id] = log_pi_tmp;
      post_params_l["log_pi_bar"] = log_pi_bar_l;
      
      for (int l = 0; l < L; l++) {
        if (node->right->id == active[l] || 
            (node->right->left  != nullptr && node->right->left->id  == active[l]) || 
            (node->right->right != nullptr && node->right->right->id == active[l])) 
        {
          active[l] = -1;
          List post_params_l = post_params[l];
          post_params_l["QTb_bar"] = NumericMatrix(n_node, d);
          post_params_l["pi_bar"] = NumericVector(n_node, 1.0 / n_node);;
          post_params_l["log_pi_bar"] = NumericVector(n_node, 0.0);
          post_params[l] = post_params_l;
        }
      }
    } else {
      // merge right child into parent
      active[right_dex] = node->id;
      merged_up[right_dex] = true;
      
      List post_params_l = post_params[right_dex];
      
      NumericMatrix QTb_bar_l = as<NumericMatrix>(post_params_l["QTb_bar"]);
      NumericVector QTb_tmp = QTb_bar_l(node->right->id, _);
      QTb_bar_l(node->right->id, _) = QTb_bar_l(node->id, _);
      QTb_bar_l(node->id, _) = QTb_tmp;
      post_params_l["QTb_bar"] = QTb_bar_l;
      
      NumericVector pi_bar_l = as<NumericVector>(post_params_l["pi_bar"]);
      double pi_tmp = pi_bar_l[node->right->id];
      pi_bar_l[node->right->id] = pi_bar_l[node->id];
      pi_bar_l[node->id] = pi_tmp;
      post_params_l["pi_bar"] = pi_bar_l;
      
      NumericVector log_pi_bar_l = as<NumericVector>(post_params_l["log_pi_bar"]);
      double log_pi_tmp = log_pi_bar_l[node->right->id];
      log_pi_bar_l[node->right->id] = log_pi_bar_l[node->id];
      log_pi_bar_l[node->id] = log_pi_tmp;
      post_params_l["log_pi_bar"] = log_pi_bar_l;
      
      post_params[right_dex] = post_params_l;
      
      if (merged_up[left_dex]) {
        // resetting left clade
        for (int l = 0; l < L; l++) {
          if (node->left->id == active[l] || 
              (node->left->left != nullptr && node->left->left->id == active[l]) || 
              (node->left->right != nullptr && node->left->right->id == active[l])) 
          {
            active[l] = -1;
            List post_params_l = post_params[l];
            post_params_l["QTb_bar"] = NumericMatrix(n_node, d);
            post_params_l["pi_bar"] = NumericVector(n_node, 1.0 / n_node);;
            post_params_l["log_pi_bar"] = NumericVector(n_node, 0.0);
            post_params[l] = post_params_l;
          }
        }
      }
    } 
  }
  return;
}

void merge_reset(
    TreeNode* root,
    NumericMatrix QTr_bar,
    List post_params,
    NumericVector lambda, 
    NumericMatrix mean_wts
) {
  int n_leaf = root->n_leaf;
  int n_node = 2 * n_leaf - 1;
  int d = lambda.length();
  int L = post_params.length();
  
  IntegerVector active(L);
  LogicalVector merged_up(L, false);
  
  List post_params_l;
  NumericVector pi_bar_l (n_node);
  NumericMatrix QTb_bar_l (n_node, d);
  NumericMatrix QTmu_bar_l (n_leaf, d);
  NumericVector muTLmu_bar_l (n_leaf, 0.0);
  
  for (int l = 0; l < L; l++) {
    double pi_max = 0.0;
    post_params_l = post_params[l];
    pi_bar_l = as<NumericVector>(post_params_l["pi_bar"]);
    for (int i=0; i < n_node; i++) {
      if (pi_bar_l[i] > pi_max) {
        pi_max = pi_bar_l[i];
        active[l] = i;
      }
    }
  }
  
  merge_reset_rec(root, n_node, d, active, merged_up, post_params); 

  // update mean params
  for (int l = 0; l < L; l++) {
    post_params_l = post_params[l];
    QTb_bar_l = as<NumericMatrix>(post_params_l["QTb_bar"]);
    QTmu_bar_l = as<NumericMatrix>(post_params_l["QTmu_bar"]);
    pi_bar_l = as<NumericVector>(post_params_l["pi_bar"]);
    
    for (int i = 0; i < n_leaf; i++) {
      for (int j = 0; j < d; j++) {
        QTr_bar(i, j) += QTmu_bar_l(i, j);
      }
    }
    
    root->mu.fill(0.0); // zero out mean before calling recursive update
    root->muTmu = 0.0;
    
    mu_bar_rec(root, lambda, mean_wts, QTmu_bar_l, muTLmu_bar_l, QTb_bar_l, pi_bar_l);
    post_params_l["QTmu_bar"] = QTmu_bar_l;
    post_params_l["muTLmu_bar"] = muTLmu_bar_l;
    post_params[l] = post_params_l;
    
    for (int i = 0; i < n_leaf; i++) {
      for (int j = 0; j < d; j++) {
        QTr_bar(i, j) -= QTmu_bar_l(i, j);
      }
    }
  }
  return;
}

NumericVector mat_tree_vb(
  NumericMatrix QTr_bar, 
  TreeNode* tree,
  int L, 
  NumericVector QTmu_0, 
  NumericVector lambda, 
  NumericMatrix Q, 
  bool fit_intercept, 
  bool fit_scale, 
  double tol, 
  double max_iter,
  bool verbose,
  NumericMatrix log_pi_l, 
  double omega_l,
  List post_params,
  NumericVector Lambda_bar_log_det, 
  NumericVector Lambda_bar_trace,
  NumericMatrix mean_wts, 
  NumericMatrix sandwich_wts
) {
  
  int n_leaf = QTr_bar.nrow();
  int d = QTr_bar.ncol();
  double d_log_omega_l = d * std::log(omega_l);
  
  List post_params_l;
  NumericMatrix QTmu_bar_l (n_leaf, d);
  
  // initialize ELBO
  std::vector<double> elbo;
  elbo.reserve(std::min(max_iter, 1e6));
  elbo.push_back(R_NegInf);
  
  // vb algorithm
  int iter = 0;
  bool scale_swich = false;
  while (iter < max_iter) {
    // update covariance matrix
    if (scale_swich){
      // todo
    }
    
    for (int l = 0; l < L; l++) {
      // add back lth partial mean
      post_params_l = post_params[l];
      QTmu_bar_l = as<NumericMatrix>(post_params_l["QTmu_bar"]);
      for (int i = 0; i < n_leaf; i++) {
        for (int j = 0; j < d; j++) {
          QTr_bar(i, j) += QTmu_bar_l(i, j);
        }
      }
      
      // update posterior parameters
      post_params_l = tree_multi_smcp(tree, lambda, QTr_bar, Lambda_bar_log_det, mean_wts, sandwich_wts, log_pi_l(_,l));
      post_params[l] = post_params_l;
      
      // subtract lth partial mean
      QTmu_bar_l = as<NumericMatrix>(post_params_l["QTmu_bar"]);
      for (int i = 0; i < n_leaf; i++) {
        for (int j = 0; j < d; j++) {
          QTr_bar(i, j) -= QTmu_bar_l(i, j);
        }
      }
    }
    
    // update intercept
    if (fit_intercept) {
      for (int j = 0; j < d; j++) {
        double mn = 0.0;
        for (int i = 0; i < n_leaf; i++) {
          QTr_bar(i, j) += QTmu_0[j];
          mn += QTr_bar(i, j);
        }
        QTmu_0[j] = mn / n_leaf;
        for (int i = 0; i < n_leaf; i++) {
          QTr_bar(i, j) -= QTmu_0[j];
        }
      }
    }
    
    // calculate elbo
    iter++;
    elbo.push_back(0.0);
    elbo[iter] = tree_multi_elbo(lambda, QTr_bar, post_params, Lambda_bar_log_det, Lambda_bar_trace, log_pi_l, omega_l, d_log_omega_l);
    
    if (verbose & (iter % 1000 == 0)) Rcout << "Iteration: " << iter << " elbo: " << elbo[iter] << "\n";
    if (std::abs((elbo[iter] - elbo[iter - 1]) / elbo[iter - 1]) < tol) break;
  }
  return wrap(elbo);
}

// [[Rcpp::export]]
List mat_mich_tree_cpp(
    NumericMatrix y, 
    IntegerMatrix edges,
    int L, 
    NumericVector mu_0, 
    NumericVector lambda, 
    NumericMatrix Q, 
    bool fit_intercept, 
    bool fit_scale, 
    double tol, 
    double max_iter,
    bool verbose,
    NumericMatrix log_pi_l, 
    double omega_l,
    List post_params
  ) {
  
  int n_leaf = y.nrow();
  int n_node = 2 * n_leaf - 1;
  int d = y.ncol();

  // construct tree
  TreeNode* tree = buildTree(edges, lambda, omega_l);
  
  List post_params_l;
  NumericMatrix QTmu_bar_l (n_leaf, d);
  NumericVector QTmu_0 (d);
  NumericVector pi_bar_l (n_node);
  
  // initialize Lambda summary params
  NumericVector Lambda_bar_log_det (n_node);
  NumericVector Lambda_bar_trace (n_node);
  NumericMatrix mean_wts (n_node, d);
  NumericMatrix sandwich_wts (n_node, d);
  
  eigen_update(tree, lambda, Lambda_bar_log_det, Lambda_bar_trace, mean_wts, sandwich_wts, omega_l);
  
  // initialize mean residual
  NumericMatrix r_bar = clone(y);
  NumericMatrix QTr_bar (n_leaf, d); // initialize decorrelated residual
  
  for (int i = 0; i < n_leaf; i++) {
    for (int j = 0; j < d; j++) {
      for (int k = 0; k < d; k++) {
        if (i == 0) QTmu_0[j] += mu_0[k] * Q(k, j);
        QTr_bar(i, j) += r_bar(i, k) * Q(k, j);
      }
      QTr_bar(i, j) -= QTmu_0[j];
    }
  }
  
  for (int l = 0; l < L; l++) {
    post_params_l = post_params[l];
    QTmu_bar_l = as<NumericMatrix>(post_params_l["QTmu_bar"]);
    for (int i = 0; i < n_leaf; i++) {
      for (int j = 0; j < d; j++) {
        QTr_bar(i, j) -= QTmu_bar_l(i, j);
      }
    }
  }
  
  // vb algorithm
  NumericVector elbo = mat_tree_vb(
    QTr_bar, tree, L, 
    QTmu_0, lambda, Q, 
    fit_intercept, fit_scale, tol, max_iter,verbose,
    log_pi_l, omega_l, 
    post_params,
    Lambda_bar_log_det, Lambda_bar_trace, mean_wts, sandwich_wts
  );
  
  merge_reset(tree, QTr_bar, post_params, lambda, mean_wts);
  
  NumericVector elbo_merged = mat_tree_vb(
    QTr_bar, tree, L, 
    QTmu_0, lambda, Q, 
    fit_intercept, fit_scale, tol, max_iter,verbose,
    log_pi_l, omega_l, 
    post_params,
    Lambda_bar_log_det, Lambda_bar_trace, mean_wts, sandwich_wts
  );
  
  bool converged = (elbo_merged.length() < max_iter);
  
  for (int i = 0; i < elbo_merged.length(); i++) {
    elbo.push_back(elbo_merged[i]);
  }
  
  // correlated intercept
  for (int i = 0; i < d; i++) {
    mu_0[i] = 0.0;
    for (int j = 0; j < d; j++) {
      mu_0[i] += QTmu_0[j] * Q(j,i);
    }
  }
  
  // creating lists of posterior parameters
  List result = List::create(
    _["L"] = L,
    _["mu_0"] = mu_0,
    _["lambda"] = lambda,
    _["Q"] = Q,   
    _["post_params"] = post_params,
    _["elbo"] = elbo,
    _["converged"] = converged
  );
  
  return result;
}