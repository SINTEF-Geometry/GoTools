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

if(WIN32)
  if(CMAKE_CL_64)
    set(WIN_LIB_DIR "win64")
#    message("The project is set to 64 bits!")
  else()
    set(WIN_LIB_DIR "win32")
#    message("The project is set to 32 bits!")
  endif()
endif()

find_library(JSONCPP_LIBRARY
  NAMES jsoncpp
#  PATHS "~/Install/jsoncpp/build/src/lib_json/Release"
  PATHS "~/Install/lib/${WIN_LIB_DIR}"
  )
