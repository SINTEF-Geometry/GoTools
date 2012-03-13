PROJECT(GoToolsCore)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoToolsCore_SOURCE_DIR}/include
  )


# Linked in libraries


# Make the gotools-core library

FILE(GLOB_RECURSE GoToolsCore_SRCS src/*.cpp src/*.C include/*.h)
ADD_LIBRARY(GoToolsCore ${GoToolsCore_SRCS})
#TARGET_LINK_LIBRARIES(GoToolsCore ${DEPLIBS})

# Apps, examples, tests, ...?
IF(GoTools_COMPILE_APPS)
  FILE(GLOB_RECURSE gotools_APPS app/*.C)
  FOREACH(app ${gotools_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoToolsCore)
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
  ENDFOREACH(app)

  FILE(GLOB_RECURSE gotools_EXAMPLES examples/*.C)
  FOREACH(app ${gotools_EXAMPLES})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoToolsCore)
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY examples)
  ENDFOREACH(app)
ENDIF(GoTools_COMPILE_APPS)

# Copy data
ADD_CUSTOM_COMMAND(
  TARGET GoToolsCore
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy_directory ${GoToolsCore_SOURCE_DIR}/data
  ${GoToolsCore_BINARY_DIR}/data
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy_directory ${GoToolsCore_SOURCE_DIR}/examples/data
  ${GoToolsCore_BINARY_DIR}/examples/data
  )


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
