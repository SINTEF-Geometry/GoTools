PROJECT(GoIgeslib)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoIgeslib_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${GoTools_COMMON_INCLUDE_DIRS}
  )


# Linked in libraries

SET(DEPLIBS
  GoToolsCore
  sisl
  )


# Make the GoIgeslib library

FILE(GLOB_RECURSE GoIgeslib_SRCS src/*.C include/*.h)
if (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoIgeslib SHARED ${GoIgeslib_SRCS})
else (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoIgeslib ${GoIgeslib_SRCS})
endif (BUILD_AS_SHARED_LIBRARY)
TARGET_LINK_LIBRARIES(GoIgeslib ${DEPLIBS})
SET_PROPERTY(TARGET GoIgeslib
  PROPERTY FOLDER "GoIgeslib/Libs")
SET_TARGET_PROPERTIES(GoIgeslib PROPERTIES SOVERSION ${GoTools_ABI_VERSION})


# Apps, examples, tests, ...?
MACRO(ADD_APPS SUBDIR PROPERTY_FOLDER IS_TEST)
  FILE(GLOB_RECURSE GoIgeslib_APPS ${SUBDIR}/*.C)
  FOREACH(app ${GoIgeslib_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoIgeslib ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${SUBDIR})
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoIgeslib/${PROPERTY_FOLDER}")
    IF(${IS_TEST})
      ADD_TEST(${appname} ${SUBDIR}/${appname}
		--log_format=XML --log_level=all --log_sink=../Testing/${appname}.xml)
      SET_TESTS_PROPERTIES( ${appname} PROPERTIES LABELS "${SUBDIR}" )
    ENDIF(${IS_TEST})
  ENDFOREACH(app)
ENDMACRO(ADD_APPS)

IF(GoTools_COMPILE_APPS)
  ADD_APPS(app "Apps" FALSE)
  ADD_APPS(examples "Examples" FALSE)
ENDIF(GoTools_COMPILE_APPS)


# Copy data
if (GoTools_COPY_DATA)
  FILE(COPY ${GoIgeslib_SOURCE_DIR}/../gotools-data/igeslib/examples/data
    DESTINATION ${GoIgeslib_BINARY_DIR}/examples)
endif()

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
  INSTALL(TARGETS GoIgeslib DESTINATION lib COMPONENT igeslib)
  # include
  INSTALL(DIRECTORY include/GoTools/igeslib
    COMPONENT igeslib-dev
    DESTINATION include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libGoIgeslib.so)
