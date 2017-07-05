#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <string>

// pll.h is missing a header guard
#ifndef LIBPLL_PLL_H_
#define LIBPLL_PLL_H_
#include <libpll/pll.h>
#endif

#include "model_parameters.hpp"
#include "pll_partition.hpp"
#include "pll_util.hpp"

TEST_CASE("partition operations are correct", "[partition]")
{
  std::string newick_path("test-data/tiny/newton.tre");
  std::string fasta_path("test-data/tiny/newton.fasta");
  std::string raxml_path("test-data/tiny/RAxML_info.newton");

  pll_utree_t* tree = pll_utree_parse_newick(newick_path.c_str());

  std::vector<std::string> labels;
  std::vector<std::string> sequences;
  pt::pll::ParseFasta(fasta_path, tree->tip_count, labels, sequences);

  pt::pll::ModelParameters parameters = pt::pll::ParseRaxmlInfo(raxml_path);

  pt::pll::Partition partition(tree, parameters, labels, sequences);

  pll_unode_t* node = tree->nodes[tree->tip_count + tree->inner_count - 1];
  partition.TraversalUpdate(node, pt::pll::TraversalType::FULL);

  SECTION("log-likelihoods are computed correctly")
  {
    REQUIRE(partition.LogLikelihood(node) == Approx(-33.387713));
  }

  SECTION("branch lengths are optimized correctly")
  {
    double original_length = node->length;
    double optimized_length = partition.OptimizeBranch(node);

    // verify that the length changed
    REQUIRE(original_length != Approx(optimized_length));

    // verify that the length is correct
    REQUIRE(optimized_length == Approx(2.607098));

    // verify that the tree was modified
    REQUIRE(node->length == optimized_length);
    REQUIRE(node->back->length == optimized_length);
  }

  SECTION("partitions initialized with a different node order are equivalent")
  {
    // clone the tree in such a way that the tips returned by
    // pll_utree_query_tipnodes() are in a different order
    pll_unode_t* other_node = pll_utree_graph_clone(node)->back->next->next;
    pll_utree_t* other_tree = pll_utree_wraptree(other_node, tree->tip_count);
    pll_utree_every(other_tree, pt::pll::cb_copy_clv_traversal);

    pt::pll::Partition other_partition(other_tree, parameters,
                                       labels, sequences);
    other_partition.TraversalUpdate(other_node, pt::pll::TraversalType::FULL);

    REQUIRE(other_partition.LogLikelihood(other_node) ==
            Approx(partition.LogLikelihood(node)));

    pll_utree_destroy(other_tree, pt::pll::cb_erase_data);
  }

  pll_utree_destroy(tree, pt::pll::cb_erase_data);
}

TEST_CASE("utility functions work correctly", "[util]")
{
  SECTION("pt::pll::ParseMultiNewick()") {
    SECTION("with a properly formed file") {
      std::string newick_path("test-data/util/two_peaks.nw");

      std::vector<pll_utree_t*> trees = pt::pll::ParseMultiNewick(newick_path);

      CHECK(trees.size() == 2);

      for (auto tree : trees) {
        REQUIRE(tree != nullptr);
        pll_utree_destroy(tree, nullptr);
      }
    }

    SECTION("with extra empty lines in file") {
      std::string newick_path("test-data/util/two_peaks_extra_line.nw");

      std::vector<pll_utree_t*> trees = pt::pll::ParseMultiNewick(newick_path);

      CHECK(trees.size() == 2);

      for (auto tree : trees) {
        REQUIRE(tree != nullptr);
        pll_utree_destroy(tree, nullptr);
      }
    }

    SECTION("with extra non-empty whitespace lines in file") {
      std::string newick_path("test-data/util/two_peaks_extra_whitespace.nw");

      REQUIRE_THROWS_AS(pt::pll::ParseMultiNewick(newick_path), std::runtime_error);
    }

    SECTION("with extra garbage in file") {
      std::string newick_path("test-data/util/two_peaks_extra_garbage.nw");

      REQUIRE_THROWS_AS(pt::pll::ParseMultiNewick(newick_path), std::runtime_error);
    }
  }
}
