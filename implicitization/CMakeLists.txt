PROJECT(GoImplicitization)

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
FIND_PACKAGE(SISL REQUIRED)
FIND_PACKAGE(Newmat REQUIRED)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoImplicitization_SOURCE_DIR}/include
  ${GoTools_INCLUDE_DIRS}
  ${SISL_INCLUDE_DIRS}
  ${Newmat_INCLUDE_DIRS}
  )


# Linked in libraries

SET(DEPLIBS
  ${SISL_LIBRARIES}
  ${Newmat_LIBRARIES}
  ${GoTools_LIBRARIES}
  )


# Make the GoImplicitization library

FILE(GLOB GoImplicitization_SRCS src/*.C)
ADD_LIBRARY(GoImplicitization ${GoImplicitization_SRCS})
TARGET_LINK_LIBRARIES(GoImplicitization ${DEPLIBS})


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

FILE(GLOB GoImplicitization_APPS app/*.C)
FILE(GLOB_RECURSE GoImplicitization_EXAMPLES examples/*.C)
FOREACH(app ${GoImplicitization_APPS} ${GoImplicitization_EXAMPLES})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} GoImplicitization ${DEPLIBS})
ENDFOREACH(app)

# Copy data
ADD_CUSTOM_COMMAND(
  TARGET GoImplicitization
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy_directory ${GoImplicitization_SOURCE_DIR}/data 
  ${GoImplicitization_BINARY_DIR}/data
  )


# 'install' target

IF(WIN32)
  # Windows
  SET(CMAKE_INSTALL_PREFIX CACHE INTERNAL "")
  SET(GoTools_INSTALL_PREFIX "$ENV{PROGRAMFILES}/SINTEF/GoTools"
    CACHE PATH "Path to install GoTools Approximate Implicitization")
  INSTALL(TARGETS GoImplicitization DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  INSTALL(DIRECTORY include/GoTools/implicitization
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS GoImplicitization DESTINATION lib/GoTools)
  # include
  INSTALL(DIRECTORY include/GoTools/implicitization
    DESTINATION include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".hpp"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)


# CPack stuff
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "GoToolsApproximateImplicitization-2.0.0")
SET(CPACK_SOURCE_IGNORE_FILES "/CVS/;/\\.svn/;\\.swp$;\\.#;/#;.*~")
INCLUDE(CPack)
