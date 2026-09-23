# SPDX-License-Identifier: MIT
find_package(Boost 1.70 CONFIG REQUIRED COMPONENTS program_options)

add_executable(legendre_convert legendre_convert.cpp)
target_link_libraries(legendre_convert PRIVATE Boost::headers Boost::program_options maxent_warnings)
install(TARGETS legendre_convert RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
