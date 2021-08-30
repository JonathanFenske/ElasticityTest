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

set(MsELA_LIBRARY_SRC
	mytools.cc
  ela_std.cc
  forces_and_parameters.cc
  postprocessing.cc)

print_all_args (
	${MsELA_LIBRARY_SRC}
)

###############################################################################
###############################################################################


###############################################################################
###############################################################################
#
# Compile and link the sources as SHARED
#

add_library (MsELA_LIBRARY SHARED ${MsELA_LIBRARY_SRC})
DEAL_II_SETUP_TARGET(MsELA_LIBRARY)


#
# Install into the DESTINATION provided by CMAKE_INSTALL_PREFIX
#
#install (TARGETS ${MsELA_LIBRARY} DESTINATION ${CMAKE_INSTALL_PREFIX})

add_executable(MsEla "main.cxx")
DEAL_II_SETUP_TARGET(MsEla)
TARGET_LINK_LIBRARIES(MsEla
	MsELA_LIBRARY)

###############################################################################
###############################################################################