PROJECT(GoTools)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)


# Instructions:
#
# 1. Comment out the lines with ADD_SUBDIRECTORY (below) you don't
# want to compile.
#
# 2. Make a build directory somewhere, and use CMake in your preferred
# fashion. Use this CMakeLists.txt as the source.
#
# 3. Compile, install, etc.
#
# jbt@sintef.no


# Do not edit this line... (When 'GoTools_BUILD_ALL' is defined, the
# Find_<Module>.cmake modules will use local include files and
# libraries instead of installed versions).
SET(GoTools_BUILD_ALL 1)


# Set various compiler flags

IF(NOT WIN32)
  # Set warning level to '-Wall' on Linux
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
ELSE(NOT WIN32)
#  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP8
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP8 /DNOMINMAX /DNOCRYPT
  -D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS")
ENDIF(NOT WIN32)

# Get rid of MESSAGEs in Release mode
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNVERBOSE")


# Extend CMake module path
OPTION(EXTEND_CMAKE_MODULE_PATH 
  "Extend the CMAKE_MODULE_PATH variable with user directories?" ON)
IF(EXTEND_CMAKE_MODULE_PATH)
  SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} 
    ${CMAKE_INSTALL_PREFIX}/share/cmake/Modules
    "./cmake/Modules"
    "~/cmake/Modules"
    "C:/cmake/Modules")
ENDIF(EXTEND_CMAKE_MODULE_PATH)

# Set output paths
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/plugins
  CACHE PATH "Path to plugin libraries")
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin
  CACHE PATH "Path to executables")

MARK_AS_ADVANCED(
  LIBRARY_OUTPUT_PATH
  EXECUTABLE_OUTPUT_PATH
  )


# Comment out the modules you don't want to build

ADD_SUBDIRECTORY(newmat)
ADD_SUBDIRECTORY(sisl)
ADD_SUBDIRECTORY(cmake)

ADD_SUBDIRECTORY(gotools-core)
ADD_SUBDIRECTORY(igeslib)
ADD_SUBDIRECTORY(parametrization)
ADD_SUBDIRECTORY(viewlib)
ADD_SUBDIRECTORY(implicitization)
ADD_SUBDIRECTORY(intersections)
ADD_SUBDIRECTORY(topology)
ADD_SUBDIRECTORY(compositemodel)
ADD_SUBDIRECTORY(qualitymodule)
ADD_SUBDIRECTORY(model_toolbox)
ADD_SUBDIRECTORY(trivariate)
ADD_SUBDIRECTORY(trivariatemodel)
ADD_SUBDIRECTORY(lr_splines)
IF(WIN32)
  ADD_SUBDIRECTORY(step_reader)
ENDIF(WIN32)

