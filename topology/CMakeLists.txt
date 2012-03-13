PROJECT(GoTopology)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoTopology_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  )


# Linked in libraries

SET(DEPLIBS
  GoToolsCore
  )


# Make the GoTopology library

FILE(GLOB_RECURSE GoTopology_SRCS src/*.C include/*.h)
ADD_LIBRARY(GoTopology ${GoTopology_SRCS})
#TARGET_LINK_LIBRARIES(GoTopology ${DEPLIBS})


# Apps, examples, tests, ...?
IF(GoTools_COMPILE_APPS)
  FILE(GLOB GoTopology_APPS app/*.C)
  FOREACH(app ${GoTopology_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoTopology ${DEPLIBS})
  ENDFOREACH(app)
ENDIF(GoTools_COMPILE_APPS)

# 'install' target

IF(WIN32)
  # Windows
  INSTALL(TARGETS GoTopology DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  INSTALL(DIRECTORY include/GoTools/topology
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  INSTALL(TARGETS GoTopology DESTINATION lib/GoTools)
  INSTALL(DIRECTORY include/GoTools/topology
    DESTINATION include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN "*.hpp"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)
