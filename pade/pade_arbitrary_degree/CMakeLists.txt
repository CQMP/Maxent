# SPDX-License-Identifier: MIT
# Built with -DMAXENT_BUILD_PADE=ON. The sources currently do not compile
# (B20 in ALPS3_MIGRATION_PLAN.md).
find_package(Boost 1.70 CONFIG REQUIRED COMPONENTS program_options)
find_path(GMP_INCLUDE_DIR gmpxx.h REQUIRED)
find_library(GMP_LIBRARY gmp REQUIRED)
find_library(GMPXX_LIBRARY gmpxx REQUIRED)

add_executable(pade pade.cpp pade_grid.cpp pade_real.cpp pade_imag.cpp pade_interpolator.cpp pade_solver.cpp)
target_include_directories(pade PRIVATE ${GMP_INCLUDE_DIR})
target_link_libraries(pade PRIVATE ${ALPSCore_LIBRARIES} Eigen3::Eigen Boost::program_options
                      ${GMPXX_LIBRARY} ${GMP_LIBRARY} maxent_warnings)
install(TARGETS pade RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
