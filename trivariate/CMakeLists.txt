PROJECT(GoTrivariate)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoTrivariate_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${sisl_SOURCE_DIR}/include
  )


# Linked in libraries

SET(DEPLIBS
  GoToolsCore
  sisl
  )


# Make the GoTrivariate library

FILE(GLOB_RECURSE GoTrivariate_SRCS src/*.C include/*.h)
ADD_LIBRARY(GoTrivariate ${GoTrivariate_SRCS})


# Apps, examples, tests, ...?

FILE(GLOB GoTrivariate_APPS app/*.C)
FOREACH(app ${GoTrivariate_APPS})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} GoTrivariate ${DEPLIBS})
  SET_TARGET_PROPERTIES(${appname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
ENDFOREACH(app)

FILE(GLOB GoTrivariate_EXAMPLES examples/*.C)
FOREACH(app ${GoTrivariate_EXAMPLES})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} GoTrivariate ${DEPLIBS})
  SET_TARGET_PROPERTIES(${appname}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY examples)
ENDFOREACH(app)


# 'install' target

IF(WIN32)
  # Windows
  # lib
  INSTALL(TARGETS GoTrivariate DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  # include
  INSTALL(DIRECTORY include/GoTools/trivariate
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS GoTrivariate DESTINATION lib/GoTools)
  # include
  INSTALL(DIRECTORY include/GoTools/trivariate
    DESTINATION include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)
