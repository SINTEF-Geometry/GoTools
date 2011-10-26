PROJECT(GoCompositeModel)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


# Find modules

# First extend CMake module path
OPTION(EXTEND_CMAKE_MODULE_PATH 
  "Extend the CMAKE_MODULE_PATH variable with user directories?" ON)
IF(EXTEND_CMAKE_MODULE_PATH)
  SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} 
    ${CMAKE_HOME_DIRECTORY}/cmake/Modules
    ${CMAKE_INSTALL_PREFIX}/share/cmake/Modules
    "./cmake/Modules"
    "~/cmake/Modules"
    "C:/cmake/Modules")
ENDIF(EXTEND_CMAKE_MODULE_PATH)

FIND_PACKAGE(Parametrization REQUIRED)
FIND_PACKAGE(GoTopology REQUIRED)
FIND_PACKAGE(GoIntersections REQUIRED)
FIND_PACKAGE(GoImplicitization REQUIRED)
FIND_PACKAGE(GoIgeslib REQUIRED)
FIND_PACKAGE(GoTools REQUIRED)
FIND_PACKAGE(TTL REQUIRED)
FIND_PACKAGE(SISL REQUIRED)
FIND_PACKAGE(Newmat REQUIRED)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoCompositeModel_SOURCE_DIR}/include
  ${Parametrization_INCLUDE_DIRS}
  ${GoTopology_INCLUDE_DIRS}
  ${GoIntersections_INCLUDE_DIRS}
  ${GoImplicitization_INCLUDE_DIRS}
  ${GoIgeslib_INCLUDE_DIRS}
  ${GoTools_INCLUDE_DIRS}
  ${TTL_INCLUDE_DIRS}
  ${SISL_INCLUDE_DIRS}
  ${Newmat_INCLUDE_DIRS}
  )


# Linked in libraries

SET(DEPLIBS
  ${Parametrization_LIBRARIES}
  ${GoTopology_LIBRARIES}
  ${GoIntersections_LIBRARIES}
  ${GoImplicitization_LIBRARIES}
  ${GoIgeslib_LIBRARIES}
  ${GoTools_LIBRARIES}
  ${TTL_LIBRARIES}
  ${SISL_LIBRARIES}
  ${Newmat_LIBRARIES}
  )


# Make the GoCompositeModel library

FILE(GLOB_RECURSE GoCompositeModel_SRCS src/*.C include/*.h)
ADD_LIBRARY(GoCompositeModel ${GoCompositeModel_SRCS})
#TARGET_LINK_LIBRARIES(GoCompositeModel ${DEPLIBS})


# Set various compiler flags

IF(NOT WIN32)
  # Set warning level to '-Wall' on Linux
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")
ELSE(NOT WIN32)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP8
  -D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS")
ENDIF(NOT WIN32)

# Get rid of MESSAGEs in Release mode
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNVERBOSE")


# Apps, examples, tests, ...?

FILE(GLOB GoCompositeModel_APPS app/*.C)
FOREACH(app ${GoCompositeModel_APPS})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} GoCompositeModel ${DEPLIBS})
ENDFOREACH(app)


## Experimental use of Gtest
#SET(GTEST_ROOT "C:/cygwin/home/jbt/prosjekter/gtest-svn")
#FIND_PATH(Gtest_INCLUDE_DIRS "gtest/gtest.h"
#  ${GTEST_ROOT}/include
#)
#INCLUDE_DIRECTORIES(
#  ${Gtest_INCLUDE_DIRS}
#)
#FIND_LIBRARY(Gtest_LIBRARIES
#  NAMES gtestd
#  PATHS ${GTEST_ROOT}/msvc/gtest-md/Debug
#)
#FILE(GLOB GoCompositeModel_TESTS tests/*.C)
#FOREACH(test ${GoCompositeModel_TESTS})
#  GET_FILENAME_COMPONENT(testname ${test} NAME_WE)
#  ADD_EXECUTABLE(${testname} ${test})
#  TARGET_LINK_LIBRARIES(${testname} GoCompositeModel ${DEPLIBS} ${Gtest_LIBRARIES})
#ENDFOREACH(test)


# Copy data
ADD_CUSTOM_COMMAND(
  TARGET GoCompositeModel
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy_directory ${GoCompositeModel_SOURCE_DIR}/data
  ${GoCompositeModel_BINARY_DIR}/data
  )


# 'install' target

IF(WIN32)
  # Windows
  SET(CMAKE_INSTALL_PREFIX CACHE INTERNAL "")
  SET(GoTools_INSTALL_PREFIX "$ENV{PROGRAMFILES}/SINTEF/GoTools"
    CACHE PATH "Path to install GoTools CompositeModel")
  # lib
  INSTALL(TARGETS GoCompositeModel DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  #include
  INSTALL(DIRECTORY include/GoTools/compositemodel
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  #lib
  INSTALL(TARGETS GoCompositeModel DESTINATION lib/GoTools)
  # include
  INSTALL(DIRECTORY include/GoTools/compositemodel
    DESTINATION include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)
