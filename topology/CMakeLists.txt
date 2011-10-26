PROJECT(GoTopology)

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

FIND_PACKAGE(GoTools REQUIRED)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoTopology_SOURCE_DIR}/include
  ${GoTools_INCLUDE_DIRS}
  )


# Linked in libraries

SET(DEPLIBS
  ${GoTools_LIBRARIES}
  )


# Make the GoTopology library

FILE(GLOB GoTopology_SRCS src/*.C)
ADD_LIBRARY(GoTopology ${GoTopology_SRCS})
TARGET_LINK_LIBRARIES(GoTopology ${DEPLIBS})


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

FILE(GLOB GoTopology_APPS app/*.C)
FOREACH(app ${GoTopology_APPS})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} GoTopology ${DEPLIBS})
ENDFOREACH(app)


# 'install' target

IF(WIN32)
  # Windows
  SET(CMAKE_INSTALL_PREFIX CACHE INTERNAL "")
  SET(GoTools_INSTALL_PREFIX "$ENV{PROGRAMFILES}/SINTEF/GoTools"
    CACHE PATH "Path to install GoTools topology")
  INSTALL(TARGETS GoTopology DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  INSTALL(DIRECTORY include/GoTools/topology
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  INSTALL(TARGETS GoTopology DESTINATION lib/GoTools)
  INSTALL(DIRECTORY include/GoTools/topology
    DESTINATION include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".hpp"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)
