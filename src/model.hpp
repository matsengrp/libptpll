#ifndef PT_PLL_MODEL_HPP_
#define PT_PLL_MODEL_HPP_

#include <string>
#include <vector>

namespace pt { namespace pll {

struct Model {
  std::string model_name;

  std::vector<double> frequencies;
  std::vector<double> subst_params;
  std::vector<double> category_rates;
};

// default number of rate categories set to 4 to match RAxML
Model ParseRaxmlInfo(const std::string& path, size_t rate_categories = 4);

} } // namespace pt::pll

#endif // PT_PLL_MODEL_HPP_
