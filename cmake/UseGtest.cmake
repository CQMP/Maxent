# Find gtest or otherwise fetch it into the build_dir/gtest
#by Andrey
if (DEFINED GTEST_ROOT)
	message(STATUS "gtest source specified at ${GTEST_ROOT}")
	find_path(GTEST_ROOT NAMES "include/gtest/gtest.h" HINTS ${GTEST_ROOT})
	if (NOT IS_DIRECTORY ${GTEST_ROOT})
            message(WARNING "Provided wrong gtest_ROOT. Please unset gtest_ROOT - gtest will be fetched")
	    unset(GTEST_ROOT CACHE)
        endif()
    endif()

    if (NOT DEFINED GTEST_ROOT)
        find_path(gtest_ROOT1 NAMES "include/gtest/gtest.h" HINTS ${CMAKE_SOURCE_DIR}/gtest-1.6.0  ${CMAKE_SOURCE_DIR}/gtest-1.7.0  ${CMAKE_SOURCE_DIR}/gtest   ${CMAKE_BINARY_DIR}/gtest)
        if (IS_DIRECTORY ${gtest_ROOT1})
		set (GTEST_ROOT ${gtest_ROOT1})
        else()
            message(STATUS "Trying to fetch gtest via subversion")
            find_package(Subversion)
            execute_process(COMMAND "${Subversion_SVN_EXECUTABLE}" "checkout" "http://googletest.googlecode.com/svn/trunk/" "gtest" WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
	    set (GTEST_ROOT "${CMAKE_BINARY_DIR}/gtest")
        endif()
    endif()

    message(STATUS "gtest is in ${GTEST_ROOT}")

    set (gtest_INCLUDE_DIR ${GTEST_ROOT}/include)
mark_as_advanced(GTEST_ROOT)
mark_as_advanced(GTEST_INCLUDE_DIR)
