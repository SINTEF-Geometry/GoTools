PROJECT(GoIgeslib)

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

FIND_PACKAGE(SISL REQUIRED)
FIND_PACKAGE(GoTools REQUIRED)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoIgeslib_SOURCE_DIR}/include
  ${GoTools_INCLUDE_DIRS}
  ${SISL_INCLUDE_DIRS}
  )


# Linked in libraries

SET(DEPLIBS
  ${GoTools_LIBRARIES}
  ${SISL_LIBRARIES}
  )


# Make the GoIgeslib library

FILE(GLOB_RECURSE GoIgeslib_SRCS src/*.C include/*.h)
ADD_LIBRARY(GoIgeslib ${GoIgeslib_SRCS})
TARGET_LINK_LIBRARIES(GoIgeslib ${DEPLIBS})


# This line may be needed only on AMD 64 bit... @jbt
# 091120: Adding this because GUV-plugins will need it. (J.O.)
IF(NOT WIN32)
  SET_TARGET_PROPERTIES(GoIgeslib PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF(NOT WIN32)


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

FILE(GLOB GoIgeslib_APPS app/*.C)
FOREACH(app ${GoIgeslib_APPS})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} GoIgeslib)
ENDFOREACH(app)


# 'install' target

IF(WIN32)
  # Windows
  SET(CMAKE_INSTALL_PREFIX CACHE INTERNAL "")
  SET(GoTools_INSTALL_PREFIX "C:/Program Files (x86)/GoTools"
    CACHE PATH "Path to install GoTools igeslib")
  # lib
  INSTALL(TARGETS GoIgeslib DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  # include
  INSTALL(DIRECTORY include/GoTools/igeslib
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS GoIgeslib DESTINATION lib/GoTools)
  # include
  INSTALL(DIRECTORY include/GoTools/igeslib
    DESTINATION include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)
