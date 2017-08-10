########################################################################
# CMake module for finding JsonCpp
#
# The following variables will be defined:
#
#  JSONCPP_FOUND
#  JSONCPP_INCLUDE_DIR
#  JSONCPP_LIBRARIES
#

find_path(JSONCPP_INCLUDE_DIR "json/json.h"
  PATHS "~/Install/jsoncpp/include"
  )
# find_library(JSONCPP_LIBRARY
#   NAMES jsoncpp libjsoncpp
#   PATHS ${LibSourcey_BUILD_DIR}/vendor/jsoncpp
# 	      "$ENV{HOME}/Install/jsoncpp/build/src/lib_json/"
#   PATH_SUFFIXES Release
#   NO_DEFAULT_PATH)


find_library(JSONCPP_LIBRARY
  NAMES jsoncpp
  PATHS "~/Install/jsoncpp/build/src/lib_json/Release"
  )

# include(${CMAKE_ROOT}/Modules/SelectLibraryConfigurations.cmake)
# select_library_configurations(JSONCPP)

# include(FindPackageHandleStandardArgs)
# find_package_handle_standard_args(JSONCPP DEFAULT_MSG JSONCPP_LIBRARIES JSONCPP_INCLUDE_DIR)
