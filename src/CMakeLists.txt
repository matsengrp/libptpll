set(PT_SRC
  model.cpp
  pll_partition.cpp
  pll_util.cpp)

add_library(ptpll_static STATIC ${PT_SRC})
add_dependencies(ptpll_static pll_modules)
target_link_libraries(ptpll_static pll_algorithm pll_optimize pll_tree pll_util pll)
