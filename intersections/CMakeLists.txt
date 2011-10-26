PROJECT(GoIntersections)

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

FIND_PACKAGE(GoImplicitization REQUIRED)
FIND_PACKAGE(GoTools REQUIRED)
FIND_PACKAGE(SISL REQUIRED)
FIND_PACKAGE(Newmat REQUIRED)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoIntersections_SOURCE_DIR}/include
  ${GoImplicitization_INCLUDE_DIRS}
  ${GoTools_INCLUDE_DIRS}
  ${SISL_INCLUDE_DIRS}
  ${Newmat_INCLUDE_DIRS}
  )


# Linked in libraries

SET(DEPLIBS
  ${GoImplicitization_LIBRARIES}
  ${GoTools_LIBRARIES}
  ${SISL_LIBRARIES}
  ${Newmat_LIBRARIES}
  )


# Make the GoIntersections library

FILE(GLOB_RECURSE GoIntersections_SRCS src/*.C include/*.h)
ADD_LIBRARY(GoIntersections ${GoIntersections_SRCS})
TARGET_LINK_LIBRARIES(GoIntersections ${DEPLIBS})


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

SET(WITH_APPS 0)
FILE(GLOB GoIntersections_APPS app/*.C)
IF(GoIntersections_APPS)
  SET(WITH_APPS 1)
ENDIF(GoIntersections_APPS)
IF(WITH_APPS)
  FIND_PACKAGE(GoIgeslib REQUIRED)
  INCLUDE_DIRECTORIES(${GoIgeslib_INCLUDE_DIRS})
  SET(DEPLIBS ${GoIgeslib_LIBRARIES} ${DEPLIBS})
ENDIF(WITH_APPS)

FOREACH(app ${GoIntersections_APPS})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} GoIntersections ${DEPLIBS})
ENDFOREACH(app)

# Copy data
IF(WITH_APPS)
  ADD_CUSTOM_COMMAND(
    TARGET GoIntersections
    POST_BUILD
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy_directory ${GoIntersections_SOURCE_DIR}/data
    ${GoIntersections_BINARY_DIR}/data
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy_directory ${GoIntersections_SOURCE_DIR}/T3_data
    ${GoIntersections_BINARY_DIR}/T3_data
    )
ENDIF(WITH_APPS)


# 'install' target

IF(WIN32)
  # Windows
  SET(CMAKE_INSTALL_PREFIX CACHE INTERNAL "")
  SET(GoTools_INSTALL_PREFIX "$ENV{PROGRAMFILES}/SINTEF/GoTools"
    CACHE PATH "Path to install GoTools Intersections")
  INSTALL(TARGETS GoIntersections DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  INSTALL(DIRECTORY include/GoTools/intersections
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS GoIntersections DESTINATION lib/GoTools)
  # include
  INSTALL(DIRECTORY include/GoTools/intersections
    DESTINATION include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".hpp"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)


# CPack stuff
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "GoToolsIntersections-2.0.0")
SET(CPACK_SOURCE_IGNORE_FILES "/CVS/;/\\.svn/;\\.swp$;\\.#;/#;.*~")
INCLUDE(CPack)
