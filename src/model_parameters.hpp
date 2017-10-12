#ifndef PT_PLL_MODEL_PARAMETERS_HPP_
#define PT_PLL_MODEL_PARAMETERS_HPP_

#include <string>
#include <vector>

namespace pt { namespace pll {

struct Model {
  std::string model_name;

  std::vector<double> frequencies;
  std::vector<double> subst_params;

  unsigned int rate_categories;
  double alpha;
};

Model ParseRaxmlInfo(const std::string& path);

} } // namespace pt::pll

#endif // PT_PLL_MODEL_PARAMETERS_HPP_
