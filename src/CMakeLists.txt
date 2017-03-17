set(PT_SRC
  model_parameters.cpp
  pll_partition.cpp
  pll_util.cpp)

add_library(ptpll SHARED ${PT_SRC})
target_link_libraries(ptpll pll)

add_library(ptpll_static STATIC ${PT_SRC})
