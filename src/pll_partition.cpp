#include "pll_partition.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <libpll/pll.h>

extern "C" {
#include <libpll/pll_optimize.h>
#include <libpll/pllmod_algorithm.h>
}

#include "model.hpp"
#include "pll_util.hpp"

namespace pt { namespace pll {

Partition::Partition(pll_utree_t* tree,
                     const Model& model,
                     const std::vector<std::string>& labels,
                     const std::vector<std::string>& sequences) :
    tip_node_count_(tree->tip_count),
    partition_(nullptr, &pll_partition_destroy),
    model_info_(nullptr, &pllmod_util_model_destroy),
    sumtable_(nullptr)
{
  // tree is already parsed

  // check that tree is healthy (i.e., every branch has a length)
  if (!TreeHealthy(tree)) {
    throw std::invalid_argument("Tree is missing branch lengths");
  }

  const unsigned int site_count = sequences[0].size();

  // TODO: replace macros/global variables with enums or arguments as
  //       appropriate. should they go in the model parameters?
  partition_ = PartitionPtr(
      pll_partition_create(tip_node_count_,     // tips
                           inner_node_count(),  // clv_buffers
                           STATES,              // states
                           site_count,          // sites
                           1,                   // rate_matrices
                           branch_count(),      // prob_matrices
                           model.category_rates.size(),  // rate_cats
                           inner_node_count(),  // scale_buffers
                           ARCH_FLAGS),         // attributes
      &pll_partition_destroy);

  SetModel(model);
  SetTipStates(tree, labels, sequences);

  params_indices_.assign(model.category_rates.size(), 0);

  AllocateScratchBuffers();
}

Partition::Partition(Partition&& rhs) :
    tip_node_count_(rhs.tip_node_count_),
    partition_(nullptr, &pll_partition_destroy),
    model_info_(nullptr, &pllmod_util_model_destroy),
    sumtable_(nullptr)
{
  partition_ = std::move(rhs.partition_);
  model_info_ = std::move(rhs.model_info_);

  params_indices_.assign(partition_->rate_cats, 0);

  AllocateScratchBuffers();
}

Partition::~Partition()
{
  FreeScratchBuffers();
}

void Partition::SetModel(const Model& model)
{
  model_info_ = ModelInfoPtr(pllmod_util_model_info_dna(model.model_name.c_str()),
                             &pllmod_util_model_destroy);

  if (!model_info_.get()) {
    throw std::invalid_argument("Invalid model name " + model.model_name);
  }

  // set frequencies at model with index 0 (we currently have only one model).
  pll_set_frequencies(partition_.get(), 0, model.frequencies.data());

  // set 6 substitution parameters at model with index 0
  pll_set_subst_params(partition_.get(), 0, model.subst_params.data());

  // set category rates
  if (model.category_rates.size() != partition_->rate_cats) {
    throw std::invalid_argument("Invalid number of category rates");
  }

  pll_set_category_rates(partition_.get(), model.category_rates.data());
}

Model Partition::GetModel() const
{
  //
  // model name
  //

  std::string model_name = model_info_->name;

  //
  // equilibrium frequencies
  //

  unsigned int states = partition_->states;

  std::vector<double> frequencies(states);
  double* p_frequencies = partition_->frequencies[0];

  std::copy(p_frequencies, p_frequencies + frequencies.size(),
            frequencies.begin());

  //
  // substitution rates
  //

  std::vector<double> subst_params((states * (states - 1)) / 2);
  double* p_subst_params = partition_->subst_params[0];

  std::copy(p_subst_params, p_subst_params + subst_params.size(),
            subst_params.begin());

  //
  // category rates
  //

  std::vector<double> category_rates(partition_->rate_cats);
  double* p_rates = partition_->rates;

  std::copy(p_rates, p_rates + category_rates.size(),
            category_rates.begin());

  return Model{model_name, frequencies, subst_params, category_rates};
}

void Partition::SetTipStates(pll_utree_t* tree,
                             const std::vector<std::string>& labels,
                             const std::vector<std::string>& sequences)
{
  if (labels.size() != sequences.size()) {
    throw std::invalid_argument("Number of labels does not match number of sequences");
  }

  if (tip_node_count_ != labels.size()) {
    throw std::invalid_argument("Unexpected number of tip nodes supplied");
  }

  std::map<std::string, unsigned int> tip_ids;

  // populate a hash table with tree tip labels
  for (unsigned int i = 0; i < tip_node_count_; ++i) {
    // the first tip_node_count_ nodes in tree->nodes are the tips
    std::string label = tree->nodes[i]->label;
    unsigned int tip_clv_index = tree->nodes[i]->clv_index;

    bool inserted;
    std::tie(std::ignore, inserted) = tip_ids.emplace(label, tip_clv_index);

    if (!inserted) {
      throw std::invalid_argument("Error inserting tip label " + label
                                  + " into map (possibly a duplicate)");
    }
  }

  // find sequences in hash table and link them with the corresponding taxa
  for (unsigned int i = 0; i < tip_node_count_; ++i) {
    auto iter = tip_ids.find(labels[i]);

    if (iter == tip_ids.end()) {
      throw std::invalid_argument("Sequence with header " + labels[i]
                                  + " does not appear in the tree");
    }

    unsigned int tip_clv_index = iter->second;
    pll_set_tip_states(partition_.get(), tip_clv_index, pll_map_nt, sequences[i].c_str());
  }
}

void Partition::AllocateScratchBuffers()
{
  // allocate scratch buffers for TraversalUpdate()
  travbuffer_.resize(node_count());
  branch_lengths_.resize(branch_count());
  matrix_indices_.resize(branch_count());
  operations_.resize(inner_node_count());

  // allocate scratch buffers for OptimizeBranch()
  sumtable_ = (double*) pll_aligned_alloc(
      partition_->sites * partition_->rate_cats * partition_->states_padded * sizeof(double),
      ALIGNMENT);
}

void Partition::FreeScratchBuffers()
{
  if (sumtable_) {
    pll_aligned_free(sumtable_);
    sumtable_ = nullptr;
  }
}

double Partition::LogLikelihood(pll_unode_t* tree, double* per_site_lnl)
{
  double lnl = pll_compute_edge_loglikelihood(
      partition_.get(), tree->clv_index, tree->scaler_index,
      tree->back->clv_index, tree->back->scaler_index, tree->pmatrix_index,
      params_indices_.data(), per_site_lnl);

  return lnl;
}

unsigned int Partition::TraversalUpdate(pll_unode_t* root, TraversalType type)
{
  unsigned int traversal_size;
  int status;

  if (type == TraversalType::FULL) {
    status = pll_utree_traverse(root, PLL_TREE_TRAVERSE_POSTORDER,
                                cb_full_traversal, travbuffer_.data(),
                                &traversal_size);
  } else if (type == TraversalType::PARTIAL) {
    status = pll_utree_traverse(root, PLL_TREE_TRAVERSE_POSTORDER,
                                cb_partial_traversal, travbuffer_.data(),
                                &traversal_size);
  } else {
    throw std::invalid_argument("Invalid traversal type");
  }

  if (!status) {
    throw std::invalid_argument("TraversalUpdate() requires an inner node");
  }

  unsigned int matrix_count;
  unsigned int ops_count;

  // Given the computed traversal descriptor, generate the operations
  // structure, and the corresponding probability matrix indices that
  // may need recomputing.
  pll_utree_create_operations(travbuffer_.data(), traversal_size,
                              branch_lengths_.data(), matrix_indices_.data(),
                              operations_.data(), &matrix_count, &ops_count);

  // Update matrix_count probability matrices for model with index 0. The i-th
  // matrix (i ranges from 0 to matrix_count - 1) is generated using branch
  // length branch_lengths[i] and can be referred to with index
  // matrix_indices[i].
  pll_update_prob_matrices(partition_.get(), params_indices_.data(),
                           matrix_indices_.data(), branch_lengths_.data(),
                           matrix_count);

  // Use the operations array to compute all ops_count inner CLVs. Operations
  // will be carried out sequentially starting from operation 0 towards
  // ops_count-1.
  pll_update_partials(partition_.get(), operations_.data(), ops_count);

  return ops_count;
}

void Partition::UpdateBranchLength(pll_unode_t* node, double length)
{
  // Update current branch lengths.
  node->length = length;
  node->back->length = length;

  // Update this branch's probability matrix now that the branch
  // length has changed. No CLVs need to be invalidated.
  pll_update_prob_matrices(partition_.get(), params_indices_.data(),
                           &(node->pmatrix_index),
                           &(node->length), 1);
}

double Partition::OptimizeBranch(pll_unode_t* node)
{
  pll_unode_t *parent = node;
  pll_unode_t *child = node->back;

  // Compute the sumtable for the particular branch once before proceeding with
  // the optimization.
  pll_update_sumtable(partition_.get(), parent->clv_index, child->clv_index,
                      parent->scaler_index, child->scaler_index,
                      params_indices_.data(), sumtable_);

  double len = node->length;
  bool maybe_decreasing = false;

  for (unsigned int i = 0; i < MAX_ITER; ++i) {
    double d1; // First derivative.
    double d2; // Second derivative.

    pll_compute_likelihood_derivatives(
        partition_.get(), parent->scaler_index, child->scaler_index, len,
        params_indices_.data(), sumtable_, &d1, &d2);

    // printf("Branch length: %f log-L: %f Derivative: %f D2: %f\n", len,
    // opt_logl, d1,d2);
    // If derivative is approximately zero then we've found the maximum.
    if (fabs(d1) < EPSILON)
      break;

    // Newton's method for finding the optimum of a function. The iteration to
    // reach the optimum is

    // x_{i+1} = x_i - f'(x_i) / f''(x_i)

    // where x_i is the current branch, f'(x_i) is the first derivative and
    // f''(x_i) is the second derivative of the likelihood function.
    if (d2 < 0.0)
      len += d1 / d2;
    else
      len -= d1 / d2;

    // If the next branch length to evaluate goes negative, we instead
    // set it to a small positive value for the next iteration. If
    // this has happened before, we stop early, as the curve is
    // probably decreasing.
    if (len < 0.0) {
      len = EPSILON;

      if (maybe_decreasing) {
        break;
      }

      maybe_decreasing = true;
    }
  }

  // No CLVs need to be invalidated; see the definition of UpdateBranchLength().
  UpdateBranchLength(node, len);

  return len;
}

void Partition::OptimizeAllBranchesOnce(pll_unode_t* tree)
{
  std::vector<pll_unode_t*> nodes(node_count(), nullptr);
  unsigned int nodes_found;

  // Traverse the entire tree and collect nodes using a callback
  // function that returns 1 for every node visited. Some of these
  // nodes will be tips, in which case we operate on node->back (the
  // tip's parent) instead of node; see below.
  if (!pll_utree_traverse(tree, PLL_TREE_TRAVERSE_POSTORDER,
                          [](pll_unode_t*) { return 1; }, nodes.data(),
                          &nodes_found)) {
    throw std::invalid_argument("OptimizeAllBranches() requires an inner node");
  }

  if (nodes_found != nodes.size()) {
    throw std::invalid_argument("Unexpected number of nodes");
  }

  for (auto node : nodes) {
    // If this is a tip node, operate on its parent instead.
    if (!node->next) {
      node = node->back;
    }

    TraversalUpdate(node, TraversalType::PARTIAL);
    OptimizeBranch(node);
  }
}

void Partition::OptimizeAllBranches(pll_unode_t* tree)
{
  TraversalUpdate(tree, TraversalType::PARTIAL);
  double loglike_prev = LogLikelihood(tree);

  OptimizeAllBranchesOnce(tree);

  TraversalUpdate(tree, TraversalType::PARTIAL);
  double loglike = LogLikelihood(tree);

  unsigned int i = 0;
  while (fabs(loglike_prev - loglike) > EPSILON && i < MAX_ITER) {
    loglike_prev = loglike;

    OptimizeAllBranchesOnce(tree);

    TraversalUpdate(tree, TraversalType::PARTIAL);
    loglike = LogLikelihood(tree);

    ++i;
  }
}

void Partition::OptimizeBranchNeighborhood(pll_unode_t* node, int radius)
{
  // TODO: passing an optimization radius of 0 to the pll-modules
  //       local optimization function appears to be broken. see
  //       https://github.com/ddarriba/pll-modules/issues/15
  if (radius <= 0) {
    throw std::invalid_argument("optimization radius must be positive");
  }

  TraversalUpdate(node, TraversalType::PARTIAL);

  // TODO: So the local branch optimization in pll-modules doesn't
  //       behave quite as I expected it to. Given the following tree:
  //
  //     a        b
  //      \      /
  //       \    /
  //      f *--- g
  //       /    \
  //      /      \
  // e___/ h      c
  //     \
  //      \
  //       d
  //
  // The starred node is the one getting passed to the optimization
  // function, and would correspond to the edge across which an NNI
  // move was made (`f->back` points at `g`). With a radius of 1, I
  // expected five branches to be optimized, `f-g`, `f-a`, `g-b`,
  // `g-c`, and `f-h`. Instead, the neighborhood appears to be defined
  // by the node itself, rather than the edge, so a radius of 1 around
  // node `f` optimizes the three edges `f-g`, `f-a`, and `f-h`. We
  // want what I expected, not what it's doing, right?
  //
  // I think I could get the expected behavior out of it by calling
  // the optimization function twice, once on `f` and once on `g`, but
  // then the edges toward the center will end up getting optimized
  // more than once, which gets to be a bigger problem as the radius
  // increases.
  //
  // For now, we're going to leave it as is, with the caveat that if
  // you want to optimize the four branches around the center,
  // you'll need to use at least a radius of 2, and accept the fact
  // that the optimization will be lopsided toward the input node.

  pllmod_opt_optimize_branch_lengths_local(partition_.get(),
                                           node,
                                           params_indices_.data(),
                                           PLLMOD_OPT_MIN_BRANCH_LEN,
                                           PLLMOD_OPT_MAX_BRANCH_LEN,
                                           EPSILON,
                                           MAX_ITER,
                                           radius,
                                           1 /* keep_update */);
}

double Partition::OptimizeModelOnce(pll_unode_t* tree)
{
  double tolerance = 1e-3;
  double lnl;

  TraversalUpdate(tree, TraversalType::PARTIAL);

  //
  // substitution rates
  //

  lnl = -1.0 * pllmod_algo_opt_subst_rates(partition_.get(),
                                           tree,
                                           0, /* params_index */
                                           params_indices_.data(),
                                           model_info_->rate_sym,
                                           PLLMOD_OPT_MIN_SUBST_RATE,
                                           PLLMOD_OPT_MAX_SUBST_RATE,
                                           PLLMOD_ALGO_BFGS_FACTR,
                                           tolerance);

  //
  // stationary frequencies
  //

  // if any of the model info's frequency symmetry values are
  // non-zero, we optimize the frequencies.
  //
  // TODO: this is a hacky substitute for
  //       pllmod_algo_opt_frequencies() not accepting the model info
  //       freq_sym array, as the pllmod_algo_opt_subst_rates()
  //       function does above with the model info rate_sym array.

  if (std::any_of(model_info_->freq_sym,
                  model_info_->freq_sym + model_info_->states,
                  [](int x) { return x != 0; }))
  {
    lnl = -1.0 * pllmod_algo_opt_frequencies(partition_.get(),
                                             tree,
                                             0, /* params_index */
                                             params_indices_.data(),
                                             PLLMOD_ALGO_BFGS_FACTR,
                                             tolerance);
  }

  //
  // gamma rate categories
  //

  if (partition_->rate_cats > 1) {
    double alpha = (PLLMOD_OPT_MIN_ALPHA + PLLMOD_OPT_MAX_ALPHA) / 2.0;

    lnl = -1.0 * pllmod_algo_opt_alpha(partition_.get(),
                                       tree,
                                       params_indices_.data(),
                                       PLLMOD_OPT_MIN_ALPHA,
                                       PLLMOD_OPT_MAX_ALPHA,
                                       &alpha,
                                       tolerance);
  }

  return lnl;
}

double Partition::OptimizeModel(pll_unode_t* tree)
{
  TraversalUpdate(tree, TraversalType::PARTIAL);
  double prev_lnl = LogLikelihood(tree);

  double curr_lnl = OptimizeModelOnce(tree);

  size_t i = 0;
  while (std::abs(prev_lnl - curr_lnl) > EPSILON && i < MAX_ITER) {
    prev_lnl = curr_lnl;
    curr_lnl = OptimizeModelOnce(tree);

    ++i;
  }

  return curr_lnl;
}

double Partition::OptimizeAllBranchesAndModel(pll_unode_t* tree)
{
  TraversalUpdate(tree, TraversalType::PARTIAL);
  double prev_lnl = LogLikelihood(tree);

  OptimizeAllBranchesOnce(tree);
  OptimizeModelOnce(tree);

  TraversalUpdate(tree, TraversalType::PARTIAL);
  double curr_lnl = LogLikelihood(tree);

  size_t i = 0;
  while (std::abs(prev_lnl - curr_lnl) > EPSILON && i < MAX_ITER) {
    prev_lnl = curr_lnl;

    OptimizeAllBranchesOnce(tree);
    OptimizeModelOnce(tree);

    // the pll-modules parameter optimization functions called in
    // OptimizeModelOnce() handle updating partials and probability
    // matrices for us, so we don't have to recompute everything with
    // a full traversal

    TraversalUpdate(tree, TraversalType::PARTIAL);
    curr_lnl = LogLikelihood(tree);

    ++i;
  }

  return curr_lnl;
}

} } // namespace pt::pll
