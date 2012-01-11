PROJECT(GoCompositeModel)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoCompositeModel_SOURCE_DIR}/include
  ${parametrization_SOURCE_DIR}/include
  ${GoTopology_SOURCE_DIR}/include
  ${GoIntersections_SOURCE_DIR}/include
  ${GoImplicitization_SOURCE_DIR}/include
  ${GoIgeslib_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${ttl_SOURCE_DIR}/include
  ${sisl_SOURCE_DIR}/include
  ${newmat_SOURCE_DIR}/include
  )


# Linked in libraries

SET(DEPLIBS
  parametrization
  GoTopology
  GoIntersections
  GoImplicitization
  GoIgeslib
  GoToolsCore
  ttl
  sisl
  newmat
  )


# Make the GoCompositeModel library

FILE(GLOB_RECURSE GoCompositeModel_SRCS src/*.C include/*.h)
ADD_LIBRARY(GoCompositeModel ${GoCompositeModel_SRCS})
#TARGET_LINK_LIBRARIES(GoCompositeModel ${DEPLIBS})


# Apps, examples, tests, ...?

FILE(GLOB_RECURSE GoCompositeModel_APPS app/*.C)
FILE(GLOB_RECURSE GoCompositeModel_EXAMPLES examples/*.C)
FOREACH(app ${GoCompositeModel_APPS} ${GoCompositeModel_EXAMPLES})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} GoCompositeModel ${DEPLIBS})
ENDFOREACH(app)


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
