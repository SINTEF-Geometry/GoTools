PROJECT(GoIsogeometricModel)

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

FIND_PACKAGE(GoTrivariateModel REQUIRED)
FIND_PACKAGE(GoCompositeModel REQUIRED)
FIND_PACKAGE(Parametrization REQUIRED)
FIND_PACKAGE(GoTopology REQUIRED)
FIND_PACKAGE(GoIntersections REQUIRED)
FIND_PACKAGE(GoImplicitization REQUIRED)
FIND_PACKAGE(GoIgeslib REQUIRED)
FIND_PACKAGE(GoTrivariate REQUIRED)
FIND_PACKAGE(GoTools REQUIRED)
FIND_PACKAGE(SISL REQUIRED)
FIND_PACKAGE(Newmat REQUIRED)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoIsogeometricModel_SOURCE_DIR}/include
  ${GoTrivariateModel_INCLUDE_DIRS}
  ${GoCompositeModel_INCLUDE_DIRS}
  ${Parametrization_INCLUDE_DIRS}
  ${GoTopology_INCLUDE_DIRS}
  ${GoIntersections_INCLUDE_DIRS}
  ${GoImplicitization_INCLUDE_DIRS}
  ${GoIgeslib_INCLUDE_DIRS}
  ${GoTrivariate_INCLUDE_DIRS}
  ${GoTools_INCLUDE_DIRS}
  ${SISL_INCLUDE_DIRS}
  ${Newmat_INCLUDE_DIRS}
  )


# Linked in libraries

SET(DEPLIBS
  ${GoTrivariateModel_LIBRARIES}
  ${GoCompositeModel_LIBRARIES}
  ${Parametrization_LIBRARIES}
  ${GoTopology_LIBRARIES}
  ${GoIntersections_LIBRARIES}
  ${GoImplicitization_LIBRARIES}
  ${GoIgeslib_LIBRARIES}
  ${GoTrivariate_LIBRARIES}
  ${GoTools_LIBRARIES}
  ${SISL_LIBRARIES}
  ${Newmat_LIBRARIES}
  )


# Make the GoIsogeometricModel library

FILE(GLOB_RECURSE GoIsogeometricModel_SRCS src/*.C include/*.h)
ADD_LIBRARY(GoIsogeometricModel ${GoIsogeometricModel_SRCS})


# Apps, examples, tests, ...?

FILE(GLOB GoIsogeometricModel_APPS app/*.C)
FOREACH(app ${GoIsogeometricModel_APPS})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} GoIsogeometricModel ${DEPLIBS})
ENDFOREACH(app)

# Copy data
ADD_CUSTOM_COMMAND(
  TARGET GoIsogeometricModel
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy_directory ${GoIsogeometricModel_SOURCE_DIR}/data 
  ${GoIsogeometricModel_BINARY_DIR}/data
  )


# 'install' target

IF(WIN32)
  # Windows
  SET(CMAKE_INSTALL_PREFIX CACHE INTERNAL "")
  SET(GoTools_INSTALL_PREFIX "$ENV{PROGRAMFILES}/SINTEF/GoTools"
    CACHE PATH "Path to install GoTools Isogeometric_model")
  # lib
  INSTALL(TARGETS GoIsogeometricModel DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  # include
  INSTALL(DIRECTORY include/GoTools/isogeometric_model
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS GoIsogeometricModel DESTINATION lib/GoTools)
  # include
  INSTALL(DIRECTORY include/GoTools/isogeometric_model
    DESTINATION include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)
