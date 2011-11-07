PROJECT(GoIntersections)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoIntersections_SOURCE_DIR}/include
  ${GoImplicitization_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${sisl_SOURCE_DIR}/include
  ${newmat_SOURCE_DIR}/include
  )


# Linked in libraries

SET(DEPLIBS
  GoImplicitization
  GoToolsCore
  sisl
  newmat
  )


# Make the GoIntersections library

FILE(GLOB_RECURSE GoIntersections_SRCS src/*.C include/*.h)
ADD_LIBRARY(GoIntersections ${GoIntersections_SRCS})
TARGET_LINK_LIBRARIES(GoIntersections ${DEPLIBS})


# Apps, examples, tests, ...?

SET(WITH_APPS 0)
FILE(GLOB GoIntersections_APPS app/*.C)
IF(GoIntersections_APPS)
  SET(WITH_APPS 1)
ENDIF(GoIntersections_APPS)
IF(WITH_APPS)
  INCLUDE_DIRECTORIES(${GoIgeslib_SOURCE_DIR}/include)
  SET(DEPLIBS GoIgeslib ${DEPLIBS})
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
    PATTERN ".hpp"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)
