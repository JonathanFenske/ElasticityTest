#############################################################
#############################################################
# Print all arguments
macro (print_all_args)
  message (STATUS "------------------------------------")
  
  set(list_var "${ARGN}")
  foreach(_currentItem IN LISTS list_var)
    message (STATUS "Adding library source file:   ${_currentItem}")
  endforeach (_currentItem)
  
  message (STATUS "------------------------------------")
endmacro (print_all_args)
#############################################################
#############################################################



###############################################################################
###############################################################################
MESSAGE(STATUS "====================================================")
MESSAGE(STATUS "============= Configuring MsELA ====================")
MESSAGE(STATUS "====================================================")

cmake_minimum_required(VERSION 2.8.12)

find_package(deal.II 9.2.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
  
if (NOT ${deal.II_FOUND})
  message(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
endif ()

IF(NOT DEAL_II_WITH_CXX11)
  MESSAGE(FATAL_ERROR "\n*** MsELA requires C++11 but your version of deal.II is not "
    "configured with it. This likely means that your compiler is too old. Check "
    "DEAL_II_WITH_CXX11 in deal.II.")
ENDIF()

IF(NOT DEAL_II_WITH_TRILINOS)
    MESSAGE(SEND_ERROR
      "\n-- deal.II was built without support for Trilinos!\n"
      )
ENDIF()

IF(NOT DEAL_II_WITH_P4EST)
    MESSAGE(SEND_ERROR
      "\n-- deal.II was built without support for p4est!\n"
      )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

message(STATUS "This is CMake ${CMAKE_VERSION} and deal.ii ${DEAL_II_VERSION}")
message(STATUS "")

###############################################################################
###############################################################################



###############################################################################
###############################################################################
ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug."
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMENT "Switch CMAKE_BUILD_TYPE to Release."
  )
  
ADD_CUSTOM_TARGET(indent
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ./contrib/utilities/indent
  COMMENT "Indenting all C++ source and header files."
  )

set(PROJECT_NAME MsELA)

# The version number
set(MsELA_VER_MAJOR 0)
set(MsELA_VER_MINOR 1)

project(${PROJECT_NAME})
###############################################################################
###############################################################################



###############################################################################
###############################################################################
# Check for the existence of various optional folders:
if (EXISTS ${CMAKE_SOURCE_DIR}/doc/CMakeLists.txt)
  set(MsELA_HAVE_DOC_DIRECTORY TRUE)
endif ()

if (EXISTS ${CMAKE_SOURCE_DIR}/test/CMakeLists.txt)
  set (MsELA_HAVE_TEST_DIRECTORY TRUE)
endif ()
###############################################################################
###############################################################################



# Change default CMAKE_INSTAL_PREFIX to ${CMAKE_BINARY_DIR}/lib
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/lib" CACHE PATH "default install path" FORCE )
endif ()

set (MsELA_LIBRARY MsELA)


###############################################################################
###############################################################################
set(MsELA_INCLUDE_DIR
		${CMAKE_SOURCE_DIR}/include)

set(MsELA_SRC_DIR 
		${CMAKE_SOURCE_DIR}/source)
add_subdirectory (${MsELA_SRC_DIR})


# Only add subdirectories for doc if exists
if (MsELA_HAVE_DOC_DIRECTORY)
  add_subdirectory (${CMAKE_SOURCE_DIR}/doc)
endif ()

# Only add subdirectories for tests if exists
if (MsELA_HAVE_TEST_DIRECTORY)
  add_subdirectory (${CMAKE_SOURCE_DIR}/test)
  ENABLE_TESTING()
endif ()
###############################################################################
###############################################################################