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

# Install the GoToolsCore target with export set
install(TARGETS parametrization
    EXPORT parametrizationTargets
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
install(EXPORT parametrizationTargets
    FILE parametrizationTargets.cmake
    #NAMESPACE GoTools::         # So consumer can do target_link_libraries(... GoTools::GoToolsCore)
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/parametrization
)

# Generate and install the Config and Version files
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/parametrizationConfigVersion.cmake"
    VERSION ${GoTools_VERSION_MAJOR}.${GoTools_VERSION_MINOR}.${GoTools_VERSION_PATCH}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/parametrizationConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/parametrizationConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/parametrization
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/parametrizationConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/parametrizationConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/parametrization
)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libparametrization.so)
