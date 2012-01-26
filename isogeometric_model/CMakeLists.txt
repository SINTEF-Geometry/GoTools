PROJECT(GoIsogeometricModel)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoIsogeometricModel_SOURCE_DIR}/include
  ${GoTrivariateModel_SOURCE_DIR}/include
  ${GoCompositeModel_SOURCE_DIR}/include
  ${parametrization_SOURCE_DIR}/include
  ${GoTopology_SOURCE_DIR}/include
  ${GoIntersections_SOURCE_DIR}/include
  ${GoImplicitization_SOURCE_DIR}/include
  ${GoIgeslib_SOURCE_DIR}/include
  ${GoTrivariate_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${sisl_SOURCE_DIR}/include
  ${newmat_SOURCE_DIR}/include
  )


# Linked in libraries

SET(DEPLIBS
  GoTrivariateModel
  GoCompositeModel
  parametrization
  GoTopology
  GoIntersections
  GoImplicitization
  GoIgeslib
  GoTrivariate
  GoToolsCore
  sisl
  newmat
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
  SET_TARGET_PROPERTIES(${appname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
ENDFOREACH(app)

FILE(GLOB GoIsogeometricModel_EXAMPLES examples/*.C)
FOREACH(app ${GoIsogeometricModel_EXAMPLES})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} GoIsogeometricModel ${DEPLIBS})
  SET_TARGET_PROPERTIES(${appname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY examples)
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
