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
TARGET_LINK_LIBRARIES(GoTopology ${DEPLIBS})
SET_PROPERTY(TARGET GoTopology
  PROPERTY FOLDER "GoTopology/Libs")


# Apps, examples, tests, ...?
IF(GoTools_COMPILE_APPS)
  FILE(GLOB GoTopology_APPS app/*.C)
  FOREACH(app ${GoTopology_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoTopology ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoTopology/Apps")
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
  INSTALL(TARGETS GoTopology DESTINATION lib COMPONENT topology)
  INSTALL(DIRECTORY include/GoTools/topology
    DESTINATION include/GoTools
    COMPONENT topology-dev
    FILES_MATCHING PATTERN "*.h"
    PATTERN "*.hpp"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libGoTopology.so)
