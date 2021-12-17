###############################################################################
###############################################################################
#
# Include directory for sources
#
include_directories(${MsELA_INCLUDE_DIR})
###############################################################################
###############################################################################


###############################################################################
###############################################################################
#
# Name all sources
#
#file(GLOB_RECURSE MsELA_TARGET_LIB_SRC  "*.cc") # source files
#set(MsELA_TARGET_LIB_SRC ${MsELA_TARGET_LIB_SRC})

set(ELA_LIBRARY_SRC
  basis_funs.cc
  lame_parameters.cc
  body_force.cc
  mytools.cc
  postprocessing.cc
  process_parameter_file.cc
  ../include/myexceptions.h)

set(MsFEM_LIBRARY_SRC
  ela_basis.cc
  ela_ms.cc
  run_problem_ms.cc)

set(StdFEM_LIBRARY_SRC
  ela_std.cc
  run_problem_std.cc)

set(COMPARE_LIBRARY_SRC
  ela_basis.cc
  ela_ms.cc
  ela_std.cc
  run_problem_compare.cc)

print_all_args (
	${ELA_LIBRARY_SRC}
  ${MsFEM_LIBRARY_SRC}
  ${StdFEM_LIBRARY_SRC}
  ${COMPARE_LIBRARY_SRC}
)

###############################################################################
###############################################################################


###############################################################################
###############################################################################
#
# Compile and link the sources as SHARED
#

add_library (ELA_LIBRARY SHARED ${ELA_LIBRARY_SRC})
DEAL_II_SETUP_TARGET(ELA_LIBRARY)

add_library(MsFEM_LIBRARY SHARED ${MsFEM_LIBRARY_SRC})
DEAL_II_SETUP_TARGET(MsFEM_LIBRARY)

add_library(StdFEM_LIBRARY SHARED ${StdFEM_LIBRARY_SRC})
DEAL_II_SETUP_TARGET(StdFEM_LIBRARY)

add_library(COMPARE_LIBRARY SHARED ${COMPARE_LIBRARY_SRC})
DEAL_II_SETUP_TARGET(COMPARE_LIBRARY)

#
# Install into the DESTINATION provided by CMAKE_INSTALL_PREFIX
#
#install (TARGETS ${MsELA_LIBRARY} DESTINATION ${CMAKE_INSTALL_PREFIX})

add_executable(MsEla "main_ms.cxx")
DEAL_II_SETUP_TARGET(MsEla)
TARGET_LINK_LIBRARIES(MsEla
  ELA_LIBRARY
  MsFEM_LIBRARY)

add_executable(StdEla "main_std.cxx")
DEAL_II_SETUP_TARGET(StdEla)
TARGET_LINK_LIBRARIES(StdEla
  ELA_LIBRARY
  StdFEM_LIBRARY)

  add_executable(compare_ela "main_compare.cxx")
  DEAL_II_SETUP_TARGET(compare_ela)
  TARGET_LINK_LIBRARIES(compare_ela
    ELA_LIBRARY
    COMPARE_LIBRARY)

###############################################################################
###############################################################################