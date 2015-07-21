# Most installations of boost will not include <boost/numeric/bindings> so this will search 
# then proceed to download it using SVN if not found
# Written by Ryan with inspiration from Andrey's script
if(DEFINED BINDINGS_ROOT)
    message(STATUS "boost/numeric/bindings specified at ${BINDINGS_ROOT}")
    find_path(BINDINGS_ROOT NAMES "boost/numeric/bindings/begin.hpp" HINTS ${BINDINGS_ROOT})
    if (NOT IS_DIRECTORY ${BINDINGS_ROOT})
	    message(WARNING "Provided wrong BINDINGS_ROOT, please remove and it will be downloaded")
	    unset(BINDINGS_ROOT CACHE)
    endif()
endif()
if(NOT DEFINED BINDINGS_ROOT)
	find_path(broot NAMES "boost/numeric/bindings/begin.hpp" HINTS ${CMAKE_SOURCE_DIR}/numeric_bindings ${CMAKE_BINARY_DIR}/numeric_bindings)
	if (IS_DIRECTORY ${broot})
		set(BINDINGS_ROOT ${broot})
	else()
		message(STATUS "Trying to fetch boost/numeric/bindings by Subversion")
		find_package(Subversion)
		execute_process(COMMAND "${Subversion_SVN_EXECUTABLE}" "checkout" "https://svn.boost.org/svn/boost/sandbox/numeric_bindings/boost/" "numeric_bindings/boost" WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
		set(BINDINGS_ROOT " ${CMAKE_BINARY_DIR}/numeric_bindings")
	endif()
endif()

message(STATUS "boost/numeric/bindings are in ${BINDINGS_ROOT}")
set(BINDINGS_INCLUDE_DIR ${BINDINGS_ROOT})

mark_as_advanced(BINDINGS_ROOT)
mark_as_advanced(BINDINGS_INCLUDE_DIR)
