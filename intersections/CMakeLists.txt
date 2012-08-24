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
  GoIgeslib
  GoImplicitization
  GoToolsCore
  sisl
  newmat
  )


# Make the GoIntersections library

FILE(GLOB_RECURSE GoIntersections_SRCS src/*.C include/*.h)
ADD_LIBRARY(GoIntersections ${GoIntersections_SRCS})
TARGET_LINK_LIBRARIES(GoIntersections ${DEPLIBS})
SET_PROPERTY(TARGET GoIntersections
  PROPERTY FOLDER "GoIntersections/Libs")


# Apps, examples, tests, ...?
IF(GoTools_COMPILE_APPS)
  FILE(GLOB GoIntersections_APPS app/*.C)
  FOREACH(app ${GoIntersections_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoIntersections ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoIntersections/Apps")
  ENDFOREACH(app)
ENDIF(GoTools_COMPILE_APPS)

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
  INSTALL(TARGETS GoIntersections DESTINATION lib COMPONENT intersections)
  # include
  INSTALL(DIRECTORY include/GoTools/intersections
    DESTINATION include/GoTools
    COMPONENT intersections-dev
    FILES_MATCHING PATTERN "*.h"
    PATTERN "*.hpp"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libGoIntersections.so)
