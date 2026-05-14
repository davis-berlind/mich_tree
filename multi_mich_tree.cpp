#include <Rcpp.h>
#include <unordered_map>
#include <vector>
using namespace Rcpp;

// Basic tree node structure
struct TreeNode {
  int id;
  double n_leaf;
  NumericMatrix::Row QTy;
  NumericVector eigen_val;
  NumericVector mean_wts;
  NumericVector sandwich_wts;
  double log_det;
  NumericMatrix Lambda;
  NumericVector mu;

  TreeNode* left = nullptr;
  TreeNode* right = nullptr;
  TreeNode* parent = nullptr;
  TreeNode* sibling = nullptr;

  TreeNode(int id_) : id(id_) {}
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

void eigen_update(TreeNode* node, NumericVector lambda, double omega_0) {
  if (node->left == nullptr && node->right == nullptr) {
    node->eigen_val = omega_0 + lambda;
  } else {
    eigen_update(node->left, lambda, omega_0);
    eigen_update(node->right, lambda, omega_0);
    node->eigen_val = node->left->eigen_val + node->right->eigen_val - omega_0; // omega_0 gets double counted
  }
  node->log_det = Rcpp::sum(Rcpp::log(node->eigen_val));
  node->mean_wts = lambda / node->eigen_val; 
  node->sandwich_wts = lambda * node->mean_wts; 
  return;
};

// Build the tree from an edge matrix
TreeNode* buildTree(
  NumericMatrix QTY,
  IntegerMatrix edges,
  NumericVector lambda,
  double omega_0
) {

  // number of leaf nodes
  int n_leaf = (QTY.nrow() + 1) / 2;
  
  // create list of nodes and children
  std::unordered_map<int, TreeNode*> nodes;
  std::unordered_map<int, std::vector<int>> children;

  for (int i = 0; i < edges.nrow(); i++) {
    // extract parent and child node index
    int parent = edges(i,0) - 1;
    int child = edges(i,1) - 1;

    children[parent].push_back(child);

    // initialize parent if necessary
    if (!nodes.count(parent)) nodes[parent] = new TreeNode(parent);
    // initialize child if necessary
    if (!nodes.count(child))  nodes[child]  = new TreeNode(child);

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

  // Find root and assign leaf data
  int root;
  for (auto& kv : nodes) {
    TreeNode* node = kv.second;
    if (node->left == nullptr && node->right == nullptr) {
      node->QTy = QTY(node->id, _);
    } else if (node->parent == nullptr) {
      root = node->id;
    }
  }

  // count leaf offspring for each node
  leaf_count(nodes[root]);
  // initialize eigen values and weights
  eigen_update(nodes[root], lambda, omega_0);

  return nodes[root];
}

void tree_multi_smcp(
  TreeNode* node,
  NumericMatrix Q,
  NumericMatrix b_bar,
  NumericMatrix QTb_bar,
  NumericVector pi_bar,
  NumericVector log_pi_bar,
  NumericVector log_pi,
  double* log_pi_max
) {
  double QTy_i, DQTy_i;
  int d = Q.ncol();

  if (node->left != nullptr || node->right != nullptr) {
    tree_smcp(node->left, Q, b_bar, QTb_bar, pi_bar, log_pi_bar, log_pi, log_pi_max);
    tree_smcp(node->right, Q, b_bar, QTb_bar, pi_bar, log_pi_bar, log_pi, log_pi_max)
    node->Qy = node->left->Qy + node->right->Qy;
  } 

  log_pi_bar[node->id] = log_pi[node->id] - 0.5 * node->log_det;
  for (int i = 0; i < d; i++) {
    QTy_i = (node->QTy)[i];
    log_pi_bar[node->id] += 0.5 * (node->sandwich_wts)[i] * QTy_i * QTy_i;
    DQTy_i = (node->mean_wts)[i] * QTy_i;
    QTb_bar(node->id, i) = DQTy_i;
    for (int j = 0; j < d; j++) {
      b_bar(node->id, j) += Q(j,i) * DQTy_i;
    }
  }

  *log_pi_max = std::max(*log_pi_max, log_pi_bar[node->id]);

  if (node->parent == nullptr) {
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

