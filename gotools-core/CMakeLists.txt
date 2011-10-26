PROJECT(GoToolsCore)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoToolsCore_SOURCE_DIR}/include
  )


# Linked in libraries

#SET(DEPLIBS ${
#  )


# Make the gotools-core library

FILE(GLOB_RECURSE GoToolsCore_SRCS src/*.cpp src/*.C include/*.h)
ADD_LIBRARY(GoToolsCore ${GoToolsCore_SRCS})
TARGET_LINK_LIBRARIES(GoToolsCore ${DEPLIBS})


# This line may be needed only on AMD 64 bit
IF(NOT WIN32)
  SET_TARGET_PROPERTIES(GoToolsCore PROPERTIES COMPILE_FLAGS "-fPIC")
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

FILE(GLOB_RECURSE gotools_APPS app/*.C)
FILE(GLOB_RECURSE gotools_EXAMPLES examples/*.C)
FOREACH(app ${gotools_APPS} ${gotools_EXAMPLES})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} GoToolsCore)
ENDFOREACH(app)


# 'install' target
IF(WIN32)
  # Windows
  SET(CMAKE_INSTALL_PREFIX CACHE INTERNAL "")
  SET(GoTools_INSTALL_PREFIX "$ENV{PROGRAMFILES}/SINTEF/GoTools"
    CACHE PATH "Path to install GoTools Core library")
  # lib
  INSTALL(TARGETS GoToolsCore DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  # include
  INSTALL(DIRECTORY include/
    DESTINATION ${GoTools_INSTALL_PREFIX}/include
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS GoToolsCore DESTINATION lib/GoTools)
  # include
  INSTALL(DIRECTORY include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)


# CPack stuff
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "GoToolsCore-2.0.1")
SET(CPACK_SOURCE_IGNORE_FILES "/CVS/;/\\.svn/;\\.swp$;\\.#;/#;.*~")
INCLUDE(CPack)
