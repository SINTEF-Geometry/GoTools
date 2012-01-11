PROJECT(GoTrivariateModel)


# Include directories

INCLUDE_DIRECTORIES(
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


# Make the GoTrivariateModel library

FILE(GLOB_RECURSE GoTrivariateModel_SRCS src/*.C include/*.h)
ADD_LIBRARY(GoTrivariateModel ${GoTrivariateModel_SRCS})
#TARGET_LINK_LIBRARIES(GoTrivariateModel ${DEPLIBS})


# Apps, examples, tests, ...?

FILE(GLOB_RECURSE GoTrivariateModel_APPS app/*.C)
FILE(GLOB_RECURSE GoTrivariateModel_EXAMPLES examples/*.C)
FOREACH(app ${GoTrivariateModel_APPS} ${GoTrivariateModel_EXAMPLES})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} GoTrivariateModel ${DEPLIBS})
ENDFOREACH(app)


# 'install' target

IF(WIN32)
  # Windows
  # lib
  INSTALL(TARGETS GoTrivariateModel DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  #include
  INSTALL(DIRECTORY include/GoTools/trivariatemodel
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  #lib
  INSTALL(TARGETS GoTrivariateModel DESTINATION lib/GoTools)
  # include
  INSTALL(DIRECTORY include/GoTools/trivariatemodel
    DESTINATION include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)
