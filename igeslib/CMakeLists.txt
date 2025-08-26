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
install(TARGETS GoIgeslib
    EXPORT GoIgeslibTargets
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
install(EXPORT GoIgeslibTargets
    FILE GoIgeslibTargets.cmake
    #NAMESPACE GoTools::         # So consumer can do target_link_libraries(... GoTools::GoToolsCore)
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GoIgeslib
)

# Generate and install the Config and Version files
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/GoIgeslibConfigVersion.cmake"
    VERSION ${GoTools_VERSION_MAJOR}.${GoTools_VERSION_MINOR}.${GoTools_VERSION_PATCH}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/GoIgeslibConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/GoIgeslibConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GoIgeslib
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/GoIgeslibConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/GoIgeslibConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GoIgeslib
)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libGoIgeslib.so)
