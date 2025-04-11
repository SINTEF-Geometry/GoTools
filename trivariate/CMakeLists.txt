PROJECT(GoTrivariate)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoTrivariate_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${GoTools_COMMON_INCLUDE_DIRS}
  )


# Linked in libraries

SET(DEPLIBS
  GoToolsCore
  sisl
  )


# Make the GoTrivariate library

FILE(GLOB_RECURSE GoTrivariate_SRCS src/*.C include/*.h)
if (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoTrivariate SHARED ${GoTrivariate_SRCS})
else (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoTrivariate ${GoTrivariate_SRCS})
endif (BUILD_AS_SHARED_LIBRARY)
TARGET_LINK_LIBRARIES(GoTrivariate ${DEPLIBS})
SET_PROPERTY(TARGET GoTrivariate
  PROPERTY FOLDER "GoTrivariate/Libs")
SET_TARGET_PROPERTIES(GoTrivariate PROPERTIES SOVERSION ${GoTools_ABI_VERSION})


# Apps, examples, tests, ...?

IF(GoTools_COMPILE_APPS)
  FILE(GLOB_RECURSE GoTrivariate_APPS app/*.C)
  FOREACH(app ${GoTrivariate_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoTrivariate ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoTrivariate/Apps")
  ENDFOREACH(app)

  FILE(GLOB_RECURSE GoTrivariate_EXAMPLES examples/*.C)
  FOREACH(app ${GoTrivariate_EXAMPLES})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoTrivariate ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY examples)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoTrivariate/Examples")
  ENDFOREACH(app)
ENDIF(GoTools_COMPILE_APPS)

# Copy data
if (GoTools_COPY_DATA)
  FILE(COPY ${GoTrivariate_SOURCE_DIR}/../gotools-data/trivariate/examples/data
    DESTINATION ${GoTrivariate_BINARY_DIR}/examples)
endif()

# === Install Targets and Export Package ===

include(CMakePackageConfigHelpers)

# 'install' target
IF(WIN32)
  # Windows
  set(CMAKE_INSTALL_LIBDIR ${GoTools_WINDOWS_INSTALL_PREFIX}/lib)
  set(CMAKE_INSTALL_INCLUDEDIR ${GoTools_WINDOWS_INSTALL_PREFIX}/include)
  set(CMAKE_INSTALL_BINDIR ${GoTools_WINDOWS_INSTALL_PREFIX}/bin)
ELSE(WIN32)
  # Linux
  include(GNUInstallDirs)
ENDIF(WIN32)

# Install the GoTrivariate target with export set
install(TARGETS GoTrivariate
    EXPORT GoTrivariateTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install headers
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h"
)

# Export the targets to a file
install(EXPORT GoTrivariateTargets
    FILE GoTrivariateTargets.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GoTrivariate
)

# Generate and install the Config and Version files
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/GoTrivariateConfigVersion.cmake"
    VERSION ${GoTools_VERSION_MAJOR}.${GoTools_VERSION_MINOR}.${GoTools_VERSION_PATCH}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/GoTrivariateConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/GoTrivariateConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GoTrivariate
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/GoTrivariateConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/GoTrivariateConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GoTrivariate
)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libGoTrivariate.so)
