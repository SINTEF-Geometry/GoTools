PROJECT(GoToolsCore)

set(CMAKE_MODULE_PATH
	${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

# OPTION(GoTools_ENABLE_OPENMP "Include OPENMP?" OFF)
# IF(GoTools_ENABLE_OPENMP)
#   FIND_PACKAGE(OpenMP REQUIRED)
# ENDIF(GoTools_ENABLE_OPENMP)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoToolsCore_SOURCE_DIR}/include
  ${GoTools_COMMON_INCLUDE_DIRS}
  )

# Linked in libraries

SET(DEPLIBS
  sisl
  )

if (GoTools_ENABLE_LOGGER)
  message("gotools-core: Defining GOTOOLS_LOG")
  find_package(spdlog REQUIRED)
  INCLUDE_DIRECTORIES(${spdlog_INCLUDE_DIR})  # Include spdlog directories
  #message("gotools-core: spdlog_INCLUDE_DIR: " ${spdlog_INCLUDE_DIR})
  IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(DEPLIBS ${DEPLIBS} ${SPDLOG_LIBRARIES} fmt)
  ENDIF()
else()
  message("Not defining GOTOOLS_LOG")
endif()

# Make the gotools-core library

FILE(GLOB_RECURSE GoToolsCore_SRCS src/*.cpp src/*.C include/*.h)
if (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoToolsCore SHARED ${GoToolsCore_SRCS})
else (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoToolsCore ${GoToolsCore_SRCS})
endif (BUILD_AS_SHARED_LIBRARY)
TARGET_LINK_LIBRARIES(GoToolsCore ${DEPLIBS})
SET_PROPERTY(TARGET GoToolsCore
  PROPERTY FOLDER "GoToolsCore/Libs")
SET_TARGET_PROPERTIES(GoToolsCore PROPERTIES SOVERSION ${GoTools_ABI_VERSION})
IF(GoTools_ENABLE_OPENMP)
  SET_TARGET_PROPERTIES(GoToolsCore PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
  SET_TARGET_PROPERTIES(GoToolsCore PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
ENDIF(GoTools_ENABLE_OPENMP)


# Apps and tests
MACRO(ADD_APPS SUBDIR PROPERTY_FOLDER IS_TEST)
  FILE(GLOB_RECURSE GoToolsCore_APPS ${SUBDIR}/*.C)
  FOREACH(app ${GoToolsCore_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoToolsCore ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${SUBDIR})
    IF(GoTools_ENABLE_OPENMP)
      SET_TARGET_PROPERTIES(${appname} PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
      SET_TARGET_PROPERTIES(${appname} PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
    ENDIF(GoTools_ENABLE_OPENMP)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoToolsCore/${PROPERTY_FOLDER}")
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

IF(GoTools_COMPILE_TESTS)
  # We check if boost-test is installed.
  SET(DEPLIBS ${DEPLIBS} ${Boost_LIBRARIES})
#  MESSAGE("Boost_LIBRARIES gotools-core: ${Boost_LIBRARIES}")
  ADD_APPS(test/unit "Unit Tests" TRUE)
  ADD_APPS(test/integration "Integration Tests" TRUE)
  ADD_APPS(test/acceptance "Acceptance Tests" TRUE)
ENDIF(GoTools_COMPILE_TESTS)


# Copy data
if (GoTools_COPY_DATA)
  FILE(COPY ${GoToolsCore_SOURCE_DIR}/../gotools-data/gotools-core/data
    DESTINATION ${GoToolsCore_BINARY_DIR})
  FILE(COPY ${GoToolsCore_SOURCE_DIR}/../gotools-data/gotools-core/examples/data
    DESTINATION ${GoToolsCore_BINARY_DIR}/examples)
endif()

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
  INSTALL(TARGETS GoToolsCore DESTINATION lib COMPONENT core)
  # include
  INSTALL(DIRECTORY include/
    DESTINATION include
    COMPONENT core-dev
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libGoToolsCore.so)
