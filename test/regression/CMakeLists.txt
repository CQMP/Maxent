# SPDX-License-Identifier: MIT
#
# Regression suite (see README.md). Two tests per set, linked by a fixture:
# labels regression-fast and regression-full (ctest -L is a regular
# expression, so -L regression selects both).
# *_generate runs the cases into the build tree, *_compare checks the results
# against reference/.
add_executable(dump_components components/dump_components.cpp)
target_link_libraries(dump_components PRIVATE maxent::core maxent_warnings)

# Prefer the python3 on PATH over macOS framework installs (e.g. Xcode's);
# override with -DPython3_EXECUTABLE=...
if(NOT DEFINED Python3_FIND_FRAMEWORK)
  set(Python3_FIND_FRAMEWORK LAST)
endif()
find_package(Python3 COMPONENTS Interpreter)
set(have_python_modules FALSE)
if(Python3_Interpreter_FOUND)
  execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy, h5py"
                  RESULT_VARIABLE python_modules_missing OUTPUT_QUIET ERROR_QUIET)
  if(python_modules_missing EQUAL 0)
    set(have_python_modules TRUE)
  endif()
endif()
if(NOT have_python_modules)
  message(WARNING "Python 3 with numpy and h5py not found: the regression tests are not registered")
  return()
endif()

add_test(NAME regression_harness_unit
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_generate.py)
set_tests_properties(regression_harness_unit PROPERTIES LABELS unit)

function(maxent_regression_test name sets label)
  set(results ${CMAKE_CURRENT_BINARY_DIR}/results-${name})
  add_test(NAME regression_${name}_generate
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/generate.py
            --maxent $<TARGET_FILE:maxent> --components $<TARGET_FILE:dump_components>
            --out ${results} --sets ${sets})
  add_test(NAME regression_${name}_compare
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare.py
            ${CMAKE_CURRENT_SOURCE_DIR}/reference ${results} --sets ${sets},components)
  set_tests_properties(regression_${name}_generate PROPERTIES
    FIXTURES_SETUP regression_${name} LABELS ${label})
  set_tests_properties(regression_${name}_compare PROPERTIES
    FIXTURES_REQUIRED regression_${name} LABELS ${label})
endfunction()

maxent_regression_test(fast fast,targeted,cli regression-fast)
if(MAXENT_BUILD_UTILITIES)
  set(results ${CMAKE_CURRENT_BINARY_DIR}/results-kk)
  add_test(NAME regression_kk_generate
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/generate.py
            --maxent $<TARGET_FILE:maxent> --kk $<TARGET_FILE:kk> --out ${results} --sets kk)
  add_test(NAME regression_kk_compare
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare.py
            ${CMAKE_CURRENT_SOURCE_DIR}/reference ${results} --sets kk)
  set_tests_properties(regression_kk_generate PROPERTIES FIXTURES_SETUP regression_kk LABELS regression-fast)
  set_tests_properties(regression_kk_compare PROPERTIES FIXTURES_REQUIRED regression_kk LABELS regression-fast)

  set(results ${CMAKE_CURRENT_BINARY_DIR}/results-legendre-convert)
  add_test(NAME regression_legendre_convert_generate
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/generate.py
            --maxent $<TARGET_FILE:maxent>
            --legendre-convert $<TARGET_FILE:legendre_convert>
            --out ${results} --sets legendre_convert)
  add_test(NAME regression_legendre_convert_compare
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare.py
            ${CMAKE_CURRENT_SOURCE_DIR}/reference ${results} --sets legendre_convert)
  set_tests_properties(regression_legendre_convert_generate PROPERTIES
    FIXTURES_SETUP regression_legendre_convert LABELS regression-fast)
  set_tests_properties(regression_legendre_convert_compare PROPERTIES
    FIXTURES_REQUIRED regression_legendre_convert LABELS regression-fast)
endif()
if(MAXENT_REGRESSION_FULL)
  maxent_regression_test(full full regression-full)
endif()
