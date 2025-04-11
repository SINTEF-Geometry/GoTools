PROJECT(GoIsogeometricModel)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoIsogeometricModel_SOURCE_DIR}/include
  ${GoTrivariateModel_SOURCE_DIR}/include
  ${GoCompositeModel_SOURCE_DIR}/include
  ${parametrization_SOURCE_DIR}/include
  ${GoTopology_SOURCE_DIR}/include
  ${GoIntersections_SOURCE_DIR}/include
  ${GoImplicitization_SOURCE_DIR}/include
  ${GoIgeslib_SOURCE_DIR}/include
  ${GoTrivariate_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${GoTools_COMMON_INCLUDE_DIRS}
  )


# Linked in libraries

SET(DEPLIBS
  GoTrivariateModel
  GoCompositeModel
  parametrization
  GoTopology
  GoIntersections
  GoImplicitization
  GoIgeslib
  GoTrivariate
  GoToolsCore
  sisl
  newmat
  )


# Make the GoIsogeometricModel library

FILE(GLOB_RECURSE GoIsogeometricModel_SRCS src/*.C include/*.h)
if (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoIsogeometricModel SHARED ${GoIsogeometricModel_SRCS})
else (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoIsogeometricModel ${GoIsogeometricModel_SRCS})
endif (BUILD_AS_SHARED_LIBRARY)
TARGET_LINK_LIBRARIES(GoIsogeometricModel ${DEPLIBS})
SET_PROPERTY(TARGET GoIsogeometricModel
  PROPERTY FOLDER "GoIsogeometricModel/Libs")
SET_TARGET_PROPERTIES(GoIsogeometricModel PROPERTIES SOVERSION ${GoTools_ABI_VERSION})


# Apps, examples, tests, ...?
IF(GoTools_COMPILE_APPS)
  FILE(GLOB GoIsogeometricModel_APPS app/*.C)
  FOREACH(app ${GoIsogeometricModel_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoIsogeometricModel ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoIsogeometricModel/Apps")
  ENDFOREACH(app)

  FILE(GLOB GoIsogeometricModel_EXAMPLES examples/*.C)
  FOREACH(app ${GoIsogeometricModel_EXAMPLES})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoIsogeometricModel ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY examples)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoIsogeometricModel/Examples")
  ENDFOREACH(app)
ENDIF(GoTools_COMPILE_APPS)

# Copy data
if (GoTools_COPY_DATA)
  ADD_CUSTOM_COMMAND(
    TARGET GoIsogeometricModel
    POST_BUILD
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy_directory ${GoIsogeometricModel_SOURCE_DIR}/../gotools-data/isogeometric_model/data
    ${GoIsogeometricModel_BINARY_DIR}/data
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
install(TARGETS GoIsogeometricModel
    EXPORT GoIsogeometricModelTargets
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
install(EXPORT GoIsogeometricModelTargets
    FILE GoIsogeometricModelTargets.cmake
    #NAMESPACE GoTools::         # So consumer can do target_link_libraries(... GoTools::GoToolsCore)
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GoIsogeometricModel
)

# Generate and install the Config and Version files
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/GoIsogeometricModelConfigVersion.cmake"
    VERSION ${GoTools_VERSION_MAJOR}.${GoTools_VERSION_MINOR}.${GoTools_VERSION_PATCH}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/GoIsogeometricModelConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/GoIsogeometricModelConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GoIsogeometricModel
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/GoIsogeometricModelConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/GoIsogeometricModelConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GoIsogeometricModel
)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libGoIsogeometricModel.so)
