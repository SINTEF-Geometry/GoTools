PROJECT(GoIntersections)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoIntersections_SOURCE_DIR}/include
  ${GoImplicitization_SOURCE_DIR}/include
  ${GoIgeslib_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${sisl_SOURCE_DIR}/include
  ${newmat_SOURCE_DIR}/include
  )


# Linked in libraries

SET(DEPLIBS
  GoImplicitization
  GoToolsCore
  GoIgeslib
  sisl
  newmat
  )


# Make the GoIntersections library

FILE(GLOB_RECURSE GoIntersections_SRCS src/*.C include/*.h)
ADD_LIBRARY(GoIntersections ${GoIntersections_SRCS})
TARGET_LINK_LIBRARIES(GoIntersections ${DEPLIBS})


# Apps, examples, tests, ...?

FILE(GLOB GoIntersections_APPS app/*.C)
FOREACH(app ${GoIntersections_APPS})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} GoIntersections ${DEPLIBS})
  SET_TARGET_PROPERTIES(${appname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
ENDFOREACH(app)


# 'install' target

IF(WIN32)
  # Windows
  # lib
  INSTALL(TARGETS GoIntersections DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  # include
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
    PATTERN "*.hpp"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)
