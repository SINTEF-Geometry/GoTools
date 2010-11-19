PROJECT(GoTools)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


# Instructions:
#
# 0. Compile and install newmat and sisl, if you haven't done so
# already.
#
# 1. Put this CMakeLists.txt file on the top-level, such that the
# source GoTools modules you would like to build are subdirectories.
#
# 2. Comment out the lines with ADD_SUBDIRECTORY (below) you don't
# want to compile.
#
# 3. Make a build directory somewhere, and use CMake in your preferred
# fashion. Use this CMakeLists.txt as the source.
#
# 4. Compile, install, etc.
#
# jbt@sintef.no


# Do not edit this line... (When 'GoTools_BUILD_ALL' is defined, the
# Find_<Module>.cmake modules will use local include files and
# libraries instead of installed versions).
SET(GoTools_BUILD_ALL 1)


# Comment out the modules you don't want to build

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
ADD_SUBDIRECTORY(step_reader)
