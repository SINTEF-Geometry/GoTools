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
  "C:/local/include"
  "/usr/include/jsoncpp"
  )

if(WIN32)
  if(${MSVC_VERSION} EQUAL 1900)
    set(MSVC_NAME "msvc2015_")
    # MESSAGE("Visual Studio 2015!")
  elseif((${MSVC_VERSION} GREATER_EQUAL 1920) AND (${MSVC_VERSION} LESS 1930))
    # MESSAGE("Visual Studio 2019!")
    set(MSVC_NAME "msvc2019_")
  elseif((${MSVC_VERSION} EQUAL 1930))
    set(MSVC_NAME "msvc2022_")
  else()
    message("MSVC version not supported or not installed!")
endif()
  if(CMAKE_CL_64)
    set(WIN_LIB_TYPE "64")
#    message("The project is set to 64 bits!")
  else()
    set(WIN_LIB_TYPE "32")
#    message("The project is set to 32 bits!")
  endif()
endif()

find_library(JSONCPP_LIBRARY_DEBUG
  NAMES jsoncpp
  PATHS "~/Install/${MSVC_NAME}lib${WIN_LIB_TYPE}/Debug"
  	"C:/local/${MSVC_NAME}lib${WIN_LIB_TYPE}/Debug"
  )

find_library(JSONCPP_LIBRARY_RELEASE
  NAMES jsoncpp
  PATHS "~/Install/${MSVC_NAME}lib${WIN_LIB_TYPE}/Release"
  	"C:/local/${MSVC_NAME}lib${WIN_LIB_TYPE}/Release"
  )

find_library(JSONCPP_LIBRARY
  NAMES jsoncpp
#  PATHS "~/Install/jsoncpp/build/src/lib_json/Release"
  PATHS "~/Install/${MSVC_NAME}lib${WIN_LIB_TYPE}"
  	"C:/local/${MSVC_NAME}lib${WIN_LIB_TYPE}"
  )

if(JSONCPP_LIBRARY_DEBUG)
  set(JSONCPP_LIBRARIES ${JSONCPP_LIBRARIES} debug ${JSONCPP_LIBRARY_DEBUG})
endif()
if(JSONCPP_LIBRARY_RELEASE)
  set(JSONCPP_LIBRARIES ${JSONCPP_LIBRARIES} optimized ${JSONCPP_LIBRARY_RELEASE})
endif()
if(JSONCPP_LIBRARY)
  set(JSONCPP_LIBRARIES ${JSONCPP_LIBRARIES} ${JSONCPP_LIBRARY})
endif()
