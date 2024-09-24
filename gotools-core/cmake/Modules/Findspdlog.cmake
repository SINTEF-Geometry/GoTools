########################################################################
# CMake module for finding spdlog
#
# The following variables will be defined:
#
#  spdlog_FOUND
#  spdlog_INCLUDE_DIR
#  spdlog_LIBRARIES
#

find_path(spdlog_INCLUDE_DIR "spdlog/spdlog.h"
#  PATHS "~/Install/spdlog/include"
  PATHS "~/Install/include"
  "C:/local/include"
  "/usr/include/spdlog"
  )

if(WIN32)
  if(${MSVC_VERSION} EQUAL 1900)
    set(MSVC_NAME "msvc2015_")
    # MESSAGE("Visual Studio 2015!")
  elseif((${MSVC_VERSION} GREATER_EQUAL 1920) AND (${MSVC_VERSION} LESS 1930))
    # MESSAGE("Visual Studio 2019!")
    set(MSVC_NAME "msvc2019_")
  elseif((${MSVC_VERSION} GREATER_EQUAL 1930) AND (${MSVC_VERSION} LESS 1950))
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

find_library(spdlog_LIBRARY_DEBUG
  NAMES spdlogd
  PATHS "~/Install/${MSVC_NAME}lib${WIN_LIB_TYPE}/Debug"
  	"C:/local/${MSVC_NAME}lib${WIN_LIB_TYPE}/Debug"
  )
#message("spdlog_LIBRARY_DEBUG: " ${spdlog_LIBRARY_DEBUG})

find_library(spdlog_LIBRARY_RELEASE
  NAMES spdlog
  PATHS "~/Install/${MSVC_NAME}lib${WIN_LIB_TYPE}/Release"
  	"C:/local/${MSVC_NAME}lib${WIN_LIB_TYPE}/Release"
  )
#message("spdlog_LIBRARY_RELEASE: " ${spdlog_LIBRARY_RELEASE})

find_library(spdlog_LIBRARY
  NAMES spdlog
  PATHS "~/Install/${MSVC_NAME}lib${WIN_LIB_TYPE}/Release"
  	"C:/local/${MSVC_NAME}lib${WIN_LIB_TYPE}/Release"
  )
#message("spdlog_LIBRARY: " ${spdlog_LIBRARY})

if(spdlog_LIBRARY_DEBUG)
  set(spdlog_LIBRARIES ${spdlog_LIBRARIES} debug ${spdlog_LIBRARY_DEBUG})
endif()
if(spdlog_LIBRARY_RELEASE)
  set(spdlog_LIBRARIES ${spdlog_LIBRARIES} optimized ${spdlog_LIBRARY_RELEASE})
endif()
if(spdlog_LIBRARY)
  set(spdlog_LIBRARIES ${spdlog_LIBRARIES} ${spdlog_LIBRARY})
endif()
#message("spdlog_LIBRARIES: " ${spdlog_LIBRARIES})
