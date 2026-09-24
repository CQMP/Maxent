# SPDX-License-Identifier: MIT
if(MAXENT_USE_SYSTEM_GTEST)
  # Prefer GTestConfig.cmake, which checks the version. Fall back to CMake's
  # FindGTest module (no version check) only if no config file exists at all.
  find_package(GTest 1.18 CONFIG QUIET)
  if(NOT GTest_FOUND)
    if(GTest_CONSIDERED_CONFIGS)
      message(FATAL_ERROR "MAXENT_USE_SYSTEM_GTEST: found GoogleTest ${GTest_CONSIDERED_VERSIONS} "
        "(${GTest_CONSIDERED_CONFIGS}), but 1.18 or newer is required.")
    endif()
    find_package(GTest REQUIRED)
    message(WARNING "MAXENT_USE_SYSTEM_GTEST: GoogleTest found without a config file; "
      "its version (1.18 or newer required) cannot be checked.")
  endif()
else()
  include(FetchContent)
  FetchContent_Declare(googletest
    URL https://github.com/google/googletest/releases/download/v1.18.0/googletest-1.18.0.tar.gz
    URL_HASH SHA256=6e3191c1455468b3fc35a417fb565c1c5071aee1b7e7f85e30cf48a98d37d8b5
    DOWNLOAD_EXTRACT_TIMESTAMP ON)
  set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endif()
include(GoogleTest)

set(unit_tests
  backcontTest
  default_modelTest
  gridTest
  legendreConvertTest
  paramFailureTest
  paramsTest
  simulationTest)

foreach(test IN LISTS unit_tests)
  add_executable(${test} ${test}.cpp)
  # GoogleTest first: the ALPSCore install ships its own (older) gtest headers in
  # its include directory, which must not shadow the ones matching this library.
  target_link_libraries(${test} PRIVATE GTest::gtest_main maxent::core maxent_warnings)
  gtest_discover_tests(${test}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    PROPERTIES LABELS unit)
endforeach()

add_subdirectory(regression)
