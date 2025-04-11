PROJECT(GoImplicitization)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoImplicitization_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${GoTools_COMMON_INCLUDE_DIRS}
  )


# Linked in libraries

SET(DEPLIBS
  GoToolsCore
  sisl
  newmat
  )


# Make the GoImplicitization library

FILE(GLOB_RECURSE GoImplicitization_SRCS src/*.C include/*.h)
if (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoImplicitization SHARED ${GoImplicitization_SRCS})
else (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoImplicitization ${GoImplicitization_SRCS})
endif (BUILD_AS_SHARED_LIBRARY)
TARGET_LINK_LIBRARIES(GoImplicitization ${DEPLIBS})
SET_PROPERTY(TARGET GoImplicitization
  PROPERTY FOLDER "GoImplicitization/Libs")
SET_TARGET_PROPERTIES(GoImplicitization PROPERTIES SOVERSION ${GoTools_ABI_VERSION})



# Apps, examples, tests, ...?
IF(GoTools_COMPILE_APPS)
  FILE(GLOB GoImplicitization_APPS app/*.C)
  FOREACH(app ${GoImplicitization_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoImplicitization ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoImplicitization/Apps")
  ENDFOREACH(app)

 FILE(GLOB_RECURSE GoImplicitization_EXAMPLES examples/*.C)
 FOREACH(app ${GoImplicitization_EXAMPLES})
   GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
   ADD_EXECUTABLE(${appname} ${app})
   TARGET_LINK_LIBRARIES(${appname} GoImplicitization ${DEPLIBS})
   SET_TARGET_PROPERTIES(${appname}
     PROPERTIES RUNTIME_OUTPUT_DIRECTORY examples)
   SET_PROPERTY(TARGET ${appname}
     PROPERTY FOLDER "GoImplicitization/Examples")
 ENDFOREACH(app)
ENDIF(GoTools_COMPILE_APPS)

# Copy data
if (GoTools_COPY_DATA)
  ADD_CUSTOM_COMMAND(
    TARGET GoImplicitization
    POST_BUILD
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy_directory ${GoImplicitization_SOURCE_DIR}/../gotools-data/implicitization/data 
    ${GoImplicitization_BINARY_DIR}/data
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
install(TARGETS GoImplicitization
    EXPORT GoImplicitizationTargets
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
install(EXPORT GoImplicitizationTargets
    FILE GoImplicitizationTargets.cmake
    #NAMESPACE GoTools::         # So consumer can do target_link_libraries(... GoTools::GoToolsCore)
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GoImplicitization
)

# Generate and install the Config and Version files
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/GoImplicitizationConfigVersion.cmake"
    VERSION ${GoTools_VERSION_MAJOR}.${GoTools_VERSION_MINOR}.${GoTools_VERSION_PATCH}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/GoImplicitizationConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/GoImplicitizationConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GoImplicitization
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/GoImplicitizationConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/GoImplicitizationConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GoImplicitization
)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libGoImplicitization.so)
