#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <string>

// pll.h is missing a header guard
#ifndef LIBPLL_PLL_H_
#define LIBPLL_PLL_H_
#include <libpll/pll.h>
#endif

#include "model.hpp"
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

  pt::pll::Model model = pt::pll::ParseRaxmlInfo(raxml_path);

  pt::pll::Partition partition(tree, model, labels, sequences);

  pll_unode_t* node = tree->nodes[tree->tip_count + tree->inner_count - 1];
  partition.TraversalUpdate(node, pt::pll::TraversalType::FULL);

  SECTION("log-likelihoods are computed correctly")
  {
    REQUIRE(partition.LogLikelihood(node) == Approx(-33.387713));
  }

  SECTION("single branch lengths are optimized correctly")
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

    pt::pll::Partition other_partition(other_tree, model,
                                       labels, sequences);
    other_partition.TraversalUpdate(other_node, pt::pll::TraversalType::FULL);

    REQUIRE(other_partition.LogLikelihood(other_node) ==
            Approx(partition.LogLikelihood(node)));

    pll_utree_destroy(other_tree, pt::pll::cb_erase_data);
  }

  pll_utree_destroy(tree, pt::pll::cb_erase_data);
}

std::map<pll_unode_t*, double> ResetBranchLengths(pll_unode_t* root,
                                                  unsigned int node_count,
                                                  double value)
{
  std::vector<pll_unode_t*> travbuffer(node_count);
  unsigned int traversal_size;

  pll_utree_traverse(root, PLL_TREE_TRAVERSE_POSTORDER,
                     [](pll_unode_t*) { return 1; }, travbuffer.data(),
                     &traversal_size);

  std::map<pll_unode_t*, double> original_lengths;

  for (auto node : travbuffer) {
    if (node != travbuffer[traversal_size - 1]->back) {
      original_lengths[node] = node->length;

      node->length = value;
      node->back->length = value;
    }
  }

  return original_lengths;
}

TEST_CASE("multiple branch lengths are optimized correctly", "[optimization]")
{
  std::string newick_path("test-data/five/RAxML_bestTree.five");
  std::string fasta_path("test-data/five/five.fasta");
  std::string raxml_path("test-data/five/RAxML_info.five");

  pll_utree_t* tree = pll_utree_parse_newick(newick_path.c_str());

  std::vector<std::string> labels;
  std::vector<std::string> sequences;
  pt::pll::ParseFasta(fasta_path, tree->tip_count, labels, sequences);

  pt::pll::Model model = pt::pll::ParseRaxmlInfo(raxml_path);

  pt::pll::Partition partition(tree, model, labels, sequences);
  pll_unode_t* root = tree->nodes[tree->tip_count + tree->inner_count - 1];

  double default_length = 1.0;
  std::map<pll_unode_t*, double> original_lengths =
      ResetBranchLengths(root, partition.node_count(), default_length);
  partition.TraversalUpdate(root, pt::pll::TraversalType::FULL);

  SECTION("branch lengths on the whole tree are optimized correctly")
  {
    partition.OptimizeAllBranches(root);

    for (auto& kv : original_lengths) {
      pll_unode_t* node = kv.first;
      double original_length = original_lengths[node];

      CHECK(node->length == Approx(original_length).epsilon(1e-3));
      CHECK(node->back->length == Approx(original_length).epsilon(1e-3));
    }
  }

  // get to a more interesting inner node
  root = root->next->back;

  SECTION("branch lengths in a neighborhood are optimized correctly") {
    SECTION("a non-positive radius is not accepted") {
      CHECK_THROWS_AS(partition.OptimizeBranchNeighborhood(root, -1),
                      std::invalid_argument);
      CHECK_THROWS_AS(partition.OptimizeBranchNeighborhood(root, 0),
                      std::invalid_argument);
    }

    SECTION("when the optimization radius is 1") {
      //std::cout << "\n\n";
      //pll_utree_show_ascii(root, PLL_UTREE_SHOW_BRANCH_LENGTH);

      partition.OptimizeBranchNeighborhood(root, 1);

      //std::cout << "\n\n";
      //pll_utree_show_ascii(root, PLL_UTREE_SHOW_BRANCH_LENGTH);

      for (auto& kv : original_lengths) {
        pll_unode_t* node = kv.first;

        /*

         |
         +---+ 0.000115
         |   |
         |   +--- 1.000000
         |   |
         |   +--- 1.000000
         |
         +--- 0.044303
         |
         +---+ 0.319134
             |
             +--- 1.000000
             |
             +--- 1.000000

        */

        double expected_length;
        if (node == root) {
          expected_length = 0.000115;
        } else if (node == root->next->back) {
          expected_length = 0.044303;
        } else if (node == root->next->next->back) {
          expected_length = 0.319134;
        } else {
          expected_length = default_length;
        }

        CHECK(node->length == Approx(expected_length));
        CHECK(node->back->length == Approx(expected_length));
      }
    }

    SECTION("when the optimization radius is 2") {
      // a radius of 2 covers the whole tree, so we expect to recover
      // all the original branch lengths

      //std::cout << "\n\n";
      //pll_utree_show_ascii(root, PLL_UTREE_SHOW_BRANCH_LENGTH);

      partition.OptimizeBranchNeighborhood(root, 2);

      //std::cout << "\n\n";
      //pll_utree_show_ascii(root, PLL_UTREE_SHOW_BRANCH_LENGTH);

      /*

       |
       +---+ 0.007016
       |   |
       |   +--- 0.053478
       |   |
       |   +--- 0.030236
       |
       +--- 0.033170
       |
       +---+ 0.074179
           |
           +--- 0.038065
           |
           +--- 0.060858

      */

      for (auto& kv : original_lengths) {
        pll_unode_t* node = kv.first;
        double original_length = original_lengths[node];

        CHECK(node->length == Approx(original_length).epsilon(1e-3));
        CHECK(node->back->length == Approx(original_length).epsilon(1e-3));
      }
    }
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

TEST_CASE("model optimization works correctly", "[model]")
{
  std::string newick_path("test-data/five/RAxML_bestTree.five");
  std::string fasta_path("test-data/five/five.fasta");
  std::string raxml_path("test-data/five/RAxML_info.five");

  pll_utree_t* tree = pll_utree_parse_newick(newick_path.c_str());

  std::vector<std::string> labels;
  std::vector<std::string> sequences;
  pt::pll::ParseFasta(fasta_path, tree->tip_count, labels, sequences);

  pt::pll::Model model = pt::pll::ParseRaxmlInfo(raxml_path);

  std::fill(model.frequencies.begin(),
            model.frequencies.end(),
            0.25);

  std::fill(model.subst_params.begin(),
            model.subst_params.end(),
            1.0);

  double alpha = 1.0;
  pll_compute_gamma_cats(alpha,
                         model.category_rates.size(),
                         model.category_rates.data(),
                         PLL_GAMMA_RATES_MEAN);

  pt::pll::Partition partition(tree, model, labels, sequences);
  pll_unode_t* root = tree->nodes[tree->tip_count + tree->inner_count - 1];
  partition.TraversalUpdate(root, pt::pll::TraversalType::FULL);

  CHECK(partition.LogLikelihood(root) == Approx(-3950.22));

  partition.OptimizeModel(root);

  CHECK(partition.LogLikelihood(root) == Approx(-3738.98));
}

TEST_CASE("full optimization works correctly", "[optimize]")
{
  std::string newick_path("test-data/five/RAxML_bestTree.five");
  std::string fasta_path("test-data/five/five.fasta");
  std::string raxml_path("test-data/five/RAxML_info.five");

  pll_utree_t* tree = pll_utree_parse_newick(newick_path.c_str());

  std::vector<std::string> labels;
  std::vector<std::string> sequences;
  pt::pll::ParseFasta(fasta_path, tree->tip_count, labels, sequences);

  pt::pll::Model model = pt::pll::ParseRaxmlInfo(raxml_path);

  std::fill(model.frequencies.begin(),
            model.frequencies.end(),
            0.25);

  std::fill(model.subst_params.begin(),
            model.subst_params.end(),
            1.0);

  double alpha = 1.0;
  pll_compute_gamma_cats(alpha,
                         model.category_rates.size(),
                         model.category_rates.data(),
                         PLL_GAMMA_RATES_MEAN);

  pt::pll::Partition partition(tree, model, labels, sequences);
  pll_unode_t* root = tree->nodes[tree->tip_count + tree->inner_count - 1];

  double default_length = 1.0;
  std::map<pll_unode_t*, double> original_lengths =
      ResetBranchLengths(root, partition.node_count(), default_length);

  partition.TraversalUpdate(root, pt::pll::TraversalType::FULL);

  CHECK(partition.LogLikelihood(root) == Approx(-6098.54));

  partition.OptimizeBranchesAndModel(root);

  CHECK(partition.LogLikelihood(root) == Approx(-3738.98));

  for (auto& kv : original_lengths) {
    pll_unode_t* node = kv.first;
    double original_length = original_lengths[node];

    CHECK(node->length == Approx(original_length).epsilon(1e-3));
    CHECK(node->back->length == Approx(original_length).epsilon(1e-3));
  }
}
