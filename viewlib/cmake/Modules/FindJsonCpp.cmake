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
#  PATHS "~/Install/jsoncpp/include"
  PATHS "~/Install/include"
  "/usr/include/jsoncpp"
  )

find_library(JSONCPP_LIBRARY
  NAMES jsoncpp
#  PATHS "~/Install/jsoncpp/build/src/lib_json/Release"
  PATHS "~/Install/lib"
  )
