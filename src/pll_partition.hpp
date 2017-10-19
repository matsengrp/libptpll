#ifndef PT_PLL_PARTITION_HPP_
#define PT_PLL_PARTITION_HPP_

#include <string>
#include <memory>
#include <vector>

#include <libpll/pll.h>

extern "C" {
#include <libpll/pllmod_util.h>
}

#include "model.hpp"
#include "pll_util.hpp"

namespace pt { namespace pll {

// MAX_ITER: Max iterations when optimizing branch lengths.
const unsigned int MAX_ITER = 32;
const double EPSILON = 1e-5; // Threshold for detecting zero.

class Partition
{
 private:
  using PartitionPtr = std::unique_ptr<pll_partition_t,
                                       decltype(&pll_partition_destroy)>;
  using ModelInfoPtr = std::unique_ptr<pllmod_subst_model_t,
                                       decltype(&pllmod_util_model_destroy)>;

  const unsigned int tip_node_count_;
  PartitionPtr partition_;
  ModelInfoPtr model_info_;
  std::vector<unsigned int> params_indices_;

  // scratch buffers for TraversalUpdate()
  std::vector<pll_unode_t*> travbuffer_;
  std::vector<double> branch_lengths_;
  std::vector<unsigned int> matrix_indices_;
  std::vector<pll_operation_t> operations_;

  // scratch buffers for OptimizeBranch()
  double* sumtable_;

 public:
  // based on old Partition constructor in partition.cpp
  Partition(pll_utree_t* tree,
            const Model& model,
            const std::vector<std::string>& labels,
            const std::vector<std::string>& sequences);

  // move constructor
  Partition(Partition&& rhs);

  // explicitly mark other operations as deleted
  Partition(const Partition& rhs) = delete;
  Partition& operator=(const Partition& rhs) = delete;
  Partition& operator=(Partition&& rhs) = delete;

  ~Partition();

  double LogLikelihood(pll_unode_t* tree, double* per_site_lnl = nullptr);

  unsigned int TraversalUpdate(pll_unode_t* root, TraversalType type);
  void UpdateBranchLength(pll_unode_t* node, double length);

  double OptimizeBranch(pll_unode_t* node);
  void OptimizeAllBranchesOnce(pll_unode_t* tree);
  void OptimizeAllBranches(pll_unode_t* tree);
  void OptimizeBranchNeighborhood(pll_unode_t* node, int radius);

  double OptimizeModelOnce(pll_unode_t* tree);
  double OptimizeModel(pll_unode_t* tree);

  double OptimizeBranchesAndModel(pll_unode_t* tree);

  Model GetModel() const;
  void SetModel(const Model& model);

  unsigned int tip_node_count() const;
  unsigned int inner_node_count() const;
  unsigned int node_count() const;
  unsigned int branch_count() const;

 private:

  // based on EquipPartitionWithData() in pll-utils.cpp
  void SetTipStates(pll_utree_t* tree, const std::vector<std::string>& labels,
                    const std::vector<std::string>& sequences);

  void AllocateScratchBuffers();
  void FreeScratchBuffers();
};

//
// Partition inlines
//

inline unsigned int Partition::tip_node_count() const
{
  return tip_node_count_;
}

inline unsigned int Partition::inner_node_count() const
{
  return tip_node_count_ - 2;
}

inline unsigned int Partition::node_count() const
{
  return tip_node_count() + inner_node_count();
}

inline unsigned int Partition::branch_count() const
{
  return node_count() - 1;
}

} } // namespace pt::pll

#endif // PT_PLL_PARTITION_HPP_
