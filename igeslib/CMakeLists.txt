PROJECT(GoIgeslib)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoIgeslib_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${sisl_SOURCE_DIR}/include
  )


# Linked in libraries

SET(DEPLIBS
  GoToolsCore
  sisl
  )


# Make the GoIgeslib library

FILE(GLOB_RECURSE GoIgeslib_SRCS src/*.C include/*.h)
ADD_LIBRARY(GoIgeslib ${GoIgeslib_SRCS})
TARGET_LINK_LIBRARIES(GoIgeslib ${DEPLIBS})


# Apps, examples, tests, ...?

FILE(GLOB GoIgeslib_APPS app/*.C)
FOREACH(app ${GoIgeslib_APPS})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} GoIgeslib)
  SET_TARGET_PROPERTIES(${appname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
ENDFOREACH(app)


# 'install' target

IF(WIN32)
  # Windows
  # lib
  INSTALL(TARGETS GoIgeslib DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  # include
  INSTALL(DIRECTORY include/GoTools/igeslib
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS GoIgeslib DESTINATION lib/GoTools)
  # include
  INSTALL(DIRECTORY include/GoTools/igeslib
    DESTINATION include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)
