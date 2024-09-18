PROJECT(parametrization)


# Include directories

INCLUDE_DIRECTORIES(
  ${parametrization_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${GoTools_COMMON_INCLUDE_DIRS}
  ${sisl_SOURCE_DIR}/include
  )


# Linked in libraries

SET(DEPLIBS
  GoToolsCore
  sisl
  )

# Make the parametrization library

FILE(GLOB_RECURSE parametrization_SRCS src/*.C include/*.h)
if (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(parametrization SHARED ${parametrization_SRCS})
else (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(parametrization ${parametrization_SRCS})
endif (BUILD_AS_SHARED_LIBRARY)
TARGET_LINK_LIBRARIES(parametrization ${DEPLIBS})
SET_PROPERTY(TARGET parametrization
  PROPERTY FOLDER "parametrization/Libs")
SET_TARGET_PROPERTIES(parametrization PROPERTIES SOVERSION ${GoTools_ABI_VERSION})


# Apps, examples, tests, ...?
IF(GoTools_COMPILE_APPS)
  FILE(GLOB_RECURSE parametrization_EXAMPLES examples/*.C)
  FOREACH(app ${parametrization_EXAMPLES})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} parametrization ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY examples)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "parametrization/Examples")
  ENDFOREACH(app)
ENDIF(GoTools_COMPILE_APPS)

# Copy data
if (GoTools_COPY_DATA)
  ADD_CUSTOM_COMMAND(
    TARGET parametrization
    POST_BUILD
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy_directory ${parametrization_SOURCE_DIR}/../gotools-data/parametrization/examples/data
    ${parametrization_BINARY_DIR}/examples/data
    )
endif()

# 'install' target
IF(WIN32)
  # Windows
  # lib
  INSTALL(TARGETS parametrization DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  # include
  INSTALL(DIRECTORY include/GoTools/parametrization
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS parametrization DESTINATION lib COMPONENT parametrization)
  # include
  INSTALL(DIRECTORY include/GoTools/parametrization
    DESTINATION include/GoTools
    COMPONENT parametrization-dev
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libparametrization.so)
