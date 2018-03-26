include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)

add_executable(test_ptpll test_ptpll.cpp)
add_dependencies(test_ptpll ptpll_static pll_modules)
target_link_libraries(test_ptpll ptpll_static lesplace-static)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test-data
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_ptpll COMMAND test_ptpll)
