#ifndef PT_PLL_UTIL_HPP_
#define PT_PLL_UTIL_HPP_

#include <string>
#include <vector>

#include <libpll/pll.h>

/// @file pll_util.hpp
/// @brief Utilities for interfacing with libpll.
/// Much of this was directly copied from libpll examples, so parts (c) Thomas
/// Flouri.

namespace pt { namespace pll {

const unsigned int STATES = 4;
const unsigned int ALIGNMENT = PLL_ALIGNMENT_AVX;
const unsigned int ARCH_FLAGS = PLL_ATTRIB_ARCH_AVX;

typedef struct { int clv_valid; } node_info_t;

enum class TraversalType { FULL, PARTIAL };

int cb_full_traversal(pll_unode_t *node);
int cb_partial_traversal(pll_unode_t *node);
int cb_copy_clv_traversal(pll_unode_t *node);
int cb_branch_healthy(pll_unode_t *tree);
void cb_erase_data(void *data);

pll_partition_t *pllext_partition_clone(pll_partition_t *rhs);

bool TreeHealthy(pll_utree_t *tree);

unsigned int ParseFasta(std::string path, unsigned int seq_count,
                        std::vector<std::string> &headers_out,
                        std::vector<std::string> &seqdata_out);

std::vector<pll_utree_t*> ParseMultiNewick(const std::string& filename);

pll_unode_t* GetVirtualRoot(pll_utree_t* tree);

void InvalidateEdgeClvs(pll_unode_t* node);

void SynchronizeTipIndices(pll_utree_t* src_tree, pll_utree_t* dest_tree);

} } // namespace pt::pll

#endif // PT_PLL_UTIL_HPP_
