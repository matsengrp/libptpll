#include "pll_util.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <libpll/pll.h>

/// @file pll_util.cpp
/// @brief Utilities for interfacing with libpll.
/// Much of this was directly copied from libpll examples, so parts (c) Thomas
/// Flouri.
/// Note that here we have a confluence of the pt C++ and the libpll C.
/// Thus this is the only place where we have mixed the two style conventions.

namespace pt { namespace pll {

// A callback function for testing if a tree has nonzero branch lengths.
int cb_branch_healthy(pll_unode_t *tree) {
  if (!tree->length)
    return 0;

  if (tree->next) {
    if (!tree->next->length || !tree->next->next->length)
      return 0;
  }
  return (tree->length == tree->back->length);
}

/// @brief Determine if a tree is healthy, i.e. has branch lengths.
/// @param[in] tree
/// A pll_utree_t to check.
/// @return Health value of tree.
bool TreeHealthy(pll_utree_t *tree) {
  return pll_utree_every(tree, cb_branch_healthy);
}

// A callback function for performing full and partial traversals.
int cb_traversal(pll_unode_t *node, TraversalType type) {
  node_info_t *node_info;

  /* if we don't want tips in the traversal we must return 0 here. For now,
     allow tips */
  if (!node->next)
    return 1;

  /* get the data element from the node and check if the CLV vector is
     oriented in the direction that we want to traverse. If the data
     element is not yet allocated then we allocate it, set the direction
     and instruct the traversal routine to place the node in the traversal array
     by returning 1 */
  node_info = (node_info_t *)(node->data);
  if (!node_info) {
    /* allocate data element */
    node->data = (node_info_t *)malloc(sizeof(node_info_t));
    node->next->data = (node_info_t *)malloc(sizeof(node_info_t));
    node->next->next->data = (node_info_t *)malloc(sizeof(node_info_t));

    /* set orientation on selected direction and traverse the subtree */
    node_info = (node_info_t *)node->data;
    node_info->clv_valid = 1;
    node_info = (node_info_t *)node->next->data;
    node_info->clv_valid = 0;
    node_info = (node_info_t *)node->next->next->data;
    node_info->clv_valid = 0;

    return 1;
  }

  /* if we want a partial traversal, the data element was already
     there, and the CLV on this direction is set, i.e. the CLV is
     valid, we instruct the traversal routine not to traverse the
     subtree rooted in this node/direction by returning 0 */
  if (node_info->clv_valid && type == TraversalType::PARTIAL)
    return 0;

  /* otherwise, set orientation on selected direction */
  node_info->clv_valid = 1;

  /* reset orientation on the other two directions and return 1 (i.e. traverse
     the subtree */
  node_info = (node_info_t *)node->next->data;
  node_info->clv_valid = 0;
  node_info = (node_info_t *)node->next->next->data;
  node_info->clv_valid = 0;

  return 1;
}

// A callback function for performing a full traversal.
int cb_full_traversal(pll_unode_t *node)
{
  return cb_traversal(node, TraversalType::FULL);
}

// A callback function for performing a partial traversal.
int cb_partial_traversal(pll_unode_t *node)
{
  return cb_traversal(node, TraversalType::PARTIAL);
}

// callback function for deep copying clv_valid values after cloning.
int cb_copy_clv_traversal(pll_unode_t *node) {

  node_info_t *node_info;

  if (!node->next)
    return 1;

  node_info = (node_info_t *)(node->data);
  if (!node_info) {
    return 1;
  }

  int node_valid = node_info->clv_valid;
  node_info = (node_info_t *)(node->next->data);
  int node_valid_1 = node_info->clv_valid;
  node_info = (node_info_t *)(node->next->next->data);
  int node_valid_2 = node_info->clv_valid;

  /* allocate data element and deep copy */
  node->data = (node_info_t *)malloc(sizeof(node_info_t));
  node_info = (node_info_t *)(node->data);
  node_info->clv_valid = node_valid;
  node->next->data = (node_info_t *)malloc(sizeof(node_info_t));
  node_info = (node_info_t *)(node->next->data);
  node_info->clv_valid = node_valid_1;
  node->next->next->data = (node_info_t *)malloc(sizeof(node_info_t));
  node_info = (node_info_t *)(node->next->next->data);
  node_info->clv_valid = node_valid_2;

  return 1;
}

// we don't use pll_utree_every() to destroy node data any more.
// pll_utree_destroy() now accepts a callback for destroying node
// data. pll_utree_destroy() will check to see if the data pointer is
// non-null as well as handle the bookkeeping necessary for freeing
// cycles of inner nodes, etc. the only thing we need to do now is
// free the data itself.
void cb_erase_data(void* data) {
    free(data);
}

/// @brief "Clone" a pll_partition_t object.
///
/// Note that this function does not make a perfect clone, but rather
/// copies the data we believe is necessary for our uses. It also does
/// not handle all the cases in which libpll can create the object.
/// For example, if ascertainment bias correction is enabled in the
/// object attributes, this function will copy CLVs incorrectly, among
/// other things.
///
/// @param[in] rhs
/// A pointer to the partition object to clone.
/// @return A pointer to the newly-cloned partition object, or nullptr on error.
pll_partition_t *pllext_partition_clone(pll_partition_t *rhs)
{
  pll_partition_t *lhs = pll_partition_create(
      rhs->tips,
      rhs->clv_buffers,
      rhs->states,
      rhs->sites,
      rhs->rate_matrices,
      rhs->prob_matrices,
      rhs->rate_cats,
      rhs->scale_buffers,
      rhs->attributes);

  if (!lhs) {
    return nullptr;
  }

  for (unsigned int i = 0; i < lhs->tips + lhs->clv_buffers; ++i) {
    for (unsigned int j = 0; j < lhs->sites * lhs->states * lhs->rate_cats; j++) {
      lhs->clv[i][j] = rhs->clv[i][j];
    }
  }

  for (unsigned int i = 0; i < lhs->prob_matrices; ++i) {
    for (unsigned int j = 0; j < lhs->states * lhs->states * lhs->rate_cats; j++) {
      lhs->pmatrix[i][j] = rhs->pmatrix[i][j];
    }
  }

  for (unsigned int j = 0; j < lhs->rate_cats; j++) {
    lhs->rates[j] = rhs->rates[j];
  }

  // TODO: this does not handle the case where rate_matrices is > 1
  // see pll.c:664
  for (unsigned int i = 0; i < lhs->states * (lhs->states - 1) / 2; i++) {
    lhs->subst_params[0][i] = rhs->subst_params[0][i];
  }

  for (unsigned int i = 0; i < lhs->rate_matrices; ++i) {
    for (unsigned int j = 0; j < lhs->states; j++) {
      lhs->frequencies[i][j] = rhs->frequencies[i][j];
    }
  }

  for (unsigned int i = 0; i < lhs->rate_matrices; ++i) {
    for (unsigned int j = 0; j < lhs->states * lhs->states; j++) {
      lhs->eigenvecs[i][j] = rhs->eigenvecs[i][j];
    }
  }

  for (unsigned int i = 0; i < lhs->rate_matrices; ++i) {
    for (unsigned int j = 0; j < lhs->states * lhs->states; j++) {
      lhs->inv_eigenvecs[i][j] = rhs->inv_eigenvecs[i][j];
    }
  }

  for (unsigned int i = 0; i < lhs->rate_matrices; ++i) {
    for (unsigned int j = 0; j < lhs->states; j++) {
      lhs->eigenvals[i][j] = rhs->eigenvals[i][j];
    }
  }

  lhs->maxstates = rhs->maxstates;

  return lhs;
}

/// @brief Parse a Fasta file.
/// @param[in] path
/// A Fasta path.
/// @param[in] seq_count
/// How many sequences are expected.
/// @param[out] headers_out
/// A vector to fill with header strings.
/// @param[out] seqdata_out
/// A vector to fill with sequence strings.
unsigned int ParseFasta(std::string path, unsigned int seq_count,
                        std::vector<std::string> &headers_out,
                        std::vector<std::string> &seqdata_out)
{
  pll_fasta_t *fp = pll_fasta_open(path.c_str(), pll_map_fasta);
  if (!fp) {
    throw std::invalid_argument("Error opening file " + path);
  }

  char *seq = NULL;
  char *hdr = NULL;
  long seqlen;
  long hdrlen;
  long seqno;

  headers_out.resize(seq_count);
  seqdata_out.resize(seq_count);

  // read FASTA sequences and make sure they are all of the same length
  unsigned int i;
  int sites = -1;
  for (i = 0; pll_fasta_getnext(fp, &hdr, &hdrlen, &seq, &seqlen, &seqno);
       ++i) {
    std::string header(hdr);
    std::string sequence(seq);
    free(hdr);
    free(seq);

    if (i >= seq_count) {
      pll_fasta_close(fp);
      throw std::invalid_argument("FASTA file contains more sequences than expected");
    }

    if (sites != -1 && sites != seqlen) {
      pll_fasta_close(fp);
      throw std::invalid_argument("FASTA file does not contain equal size sequences");
    }

    if (sites == -1)
      sites = seqlen;

    headers_out[i] = header;
    seqdata_out[i] = sequence;
  }

  // close FASTA file
  pll_fasta_close(fp);

  // did we stop reading the file because we reached EOF?
  if (pll_errno != PLL_ERROR_FILE_EOF) {
    throw std::runtime_error("Error while reading file " + path);
  }

  if (sites < 0) {
    throw std::runtime_error("Unable to read alignment");
  }

  if (i != seq_count) {
    throw std::invalid_argument("Some taxa are missing from FASTA file");
  }

  return (unsigned int)sites;
}

std::vector<pll_utree_t*> ParseMultiNewick(const std::string& filename)
{
  std::ifstream ifs(filename);

  if (!ifs) {
    throw std::invalid_argument("Error opening file " + filename);
  }

  // this function can handle empty lines in the multi-Newick file,
  // but anything else that fails to be parsed by libpll will result
  // in an exception being thrown.

  std::vector<pll_utree_t*> trees;
  for (std::string line; std::getline(ifs, line); ) {
    if (line.empty()) {
      continue;
    }

    pll_utree_t* tree = pll_utree_parse_newick_string(line.c_str());

    if (!tree) {
      for (auto t : trees) {
        pll_utree_destroy(t, nullptr);
      }

      throw std::runtime_error("Failed to parse '" + line +
                               "' as a Newick string");
    }

    trees.emplace_back(tree);
  }

  // partitions maintain state related to the tree tip node/partition
  // data indices, which means that if we want to use the same
  // partition for evaluating these trees, their tips need to be
  // synchronized based on the tip labels.

  for (size_t i = 1; i < trees.size(); ++i) {
    SynchronizeTipIndices(trees[0], trees[i]);
  }

  return trees;
}

pll_unode_t* GetVirtualRoot(pll_utree_t* tree)
{
  return tree->nodes[tree->tip_count + tree->inner_count - 1];
}

void InvalidateEdgeClvs(pll_unode_t* node)
{
  node_info_t* node_info;

  if (node->data) {
    node_info = static_cast<node_info_t*>(node->data);
    node_info->clv_valid = 0;
  }

  if (node->back->data) {
    node_info = static_cast<node_info_t*>(node->back->data);
    node_info->clv_valid = 0;
  }
}

void SynchronizeTipIndices(pll_utree_t* src_tree, pll_utree_t* dest_tree)
{
  if (src_tree->tip_count != dest_tree->tip_count) {
    throw std::invalid_argument("Number of tips differs in source and destination");
  }

  std::map<std::string, pll_unode_t*> src_tips;

  for (unsigned int i = 0; i < src_tree->tip_count; ++i) {
    pll_unode_t* node = src_tree->nodes[i];
    std::string label = node->label;

    bool inserted;
    std::tie(std::ignore, inserted) = src_tips.emplace(label, node);

    if (!inserted) {
      throw std::invalid_argument("Error inserting tip label " + label
                                  + " into map (possibly a duplicate)");
    }
  }

  for (unsigned int i = 0; i < src_tree->tip_count; ++i) {
    pll_unode_t* dest_node = dest_tree->nodes[i];
    std::string label = dest_node->label;

    auto iter = src_tips.find(label);
    if (iter == src_tips.end()) {
      throw std::invalid_argument("Tip with label " + label +
                                  " does not appear in source tree");
    }

    pll_unode_t* src_node = iter->second;

    dest_node->node_index = src_node->node_index;
    dest_node->clv_index = src_node->clv_index;
    dest_node->pmatrix_index = src_node->pmatrix_index;
    dest_node->scaler_index = src_node->scaler_index;
  }
}

} } // namespace pt::pll
