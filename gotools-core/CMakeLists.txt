PROJECT(GoToolsCore)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoToolsCore_SOURCE_DIR}/include
  )


# Linked in libraries

SET(DEPLIBS
  )


# Make the gotools-core library

FILE(GLOB_RECURSE GoToolsCore_SRCS src/*.cpp src/*.C include/*.h)
ADD_LIBRARY(GoToolsCore ${GoToolsCore_SRCS})
#TARGET_LINK_LIBRARIES(GoToolsCore ${DEPLIBS})

# This line may be needed only on AMD 64 bit
IF(NOT WIN32)
  SET_TARGET_PROPERTIES(GoToolsCore PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF(NOT WIN32)


# Apps, examples, tests, ...?

FILE(GLOB_RECURSE gotools_APPS app/*.C)
FILE(GLOB_RECURSE gotools_EXAMPLES examples/*.C)
FOREACH(app ${gotools_APPS} ${gotools_EXAMPLES})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} GoToolsCore)
ENDFOREACH(app)


# 'install' target
IF(WIN32)
  # Windows
  # lib
  INSTALL(TARGETS GoToolsCore DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  # include
  INSTALL(DIRECTORY include/
    DESTINATION ${GoTools_INSTALL_PREFIX}/include
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS GoToolsCore DESTINATION lib/GoTools)
  # include
  INSTALL(DIRECTORY include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)
