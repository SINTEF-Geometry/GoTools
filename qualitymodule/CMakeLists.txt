PROJECT(GoQualityModule)

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

FIND_PACKAGE(GoCompositeModel REQUIRED)
FIND_PACKAGE(Parametrization REQUIRED)
FIND_PACKAGE(GoTopology REQUIRED)
FIND_PACKAGE(GoIntersections REQUIRED)
FIND_PACKAGE(GoImplicitization REQUIRED)
FIND_PACKAGE(GoIgeslib REQUIRED)
FIND_PACKAGE(GoTools REQUIRED)
FIND_PACKAGE(SISL REQUIRED)
FIND_PACKAGE(Newmat REQUIRED)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoQualityModule_SOURCE_DIR}/include
  ${GoCompositeModel_INCLUDE_DIRS}
  ${Parametrization_INCLUDE_DIRS}
  ${GoTopology_INCLUDE_DIRS}
  ${GoIntersections_INCLUDE_DIRS}
  ${GoImplicitization_INCLUDE_DIRS}
  ${GoIgeslib_INCLUDE_DIRS}
  ${GoTools_INCLUDE_DIRS}
  ${SISL_INCLUDE_DIRS}
  ${Newmat_INCLUDE_DIRS}
  )


# Linked in libraries

SET(DEPLIBS
  ${GoCompositeModel_LIBRARIES}
  ${Parametrization_LIBRARIES}
  ${GoTopology_LIBRARIES}
  ${GoIntersections_LIBRARIES}
  ${GoImplicitization_LIBRARIES}
  ${GoIgeslib_LIBRARIES}
  ${GoTools_LIBRARIES}
  ${SISL_LIBRARIES}
  ${Newmat_LIBRARIES}
  )


# Make the GoQualityModule library

FILE(GLOB_RECURSE GoQualityModule_SRCS src/*.C include/*.h)
ADD_LIBRARY(GoQualityModule ${GoQualityModule_SRCS})
#TARGET_LINK_LIBRARIES(GoQualityModule ${DEPLIBS})


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

FILE(GLOB GoQualityModule_APPS app/*.C)
FOREACH(app ${GoQualityModule_APPS})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} GoQualityModule ${DEPLIBS})
ENDFOREACH(app)

# Copy data
ADD_CUSTOM_COMMAND(
  TARGET GoQualityModule
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy_directory ${GoQualityModule_SOURCE_DIR}/data
  ${GoQualityModule_BINARY_DIR}/data
  )


# 'install' target

IF(WIN32)
  # Windows
  SET(CMAKE_INSTALL_PREFIX CACHE INTERNAL "")
  SET(GoTools_INSTALL_PREFIX "$ENV{PROGRAMFILES}/SINTEF/GoTools"
    CACHE PATH "Path to install GoTools QualityModule")
  # lib
  INSTALL(TARGETS GoQualityModule DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  # include
  INSTALL(DIRECTORY include/GoTools/qualitymodule
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS GoQualityModule DESTINATION lib/GoTools)
  # include
  INSTALL(DIRECTORY include/GoTools/qualitymodule
    DESTINATION include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)
