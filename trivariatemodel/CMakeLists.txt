PROJECT(GoTrivariateModel)

set(CMAKE_MODULE_PATH
	${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

# Find modules

FIND_PACKAGE(PugiXML REQUIRED)

# Include directories

INCLUDE_DIRECTORIES(
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
  ${PUGIXML_INCLUDE_DIR}
  )


# Linked in libraries

SET(DEPLIBS
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
  ${PUGIXML_LIBRARIES}
  )


# Make the GoTrivariateModel library

FILE(GLOB_RECURSE GoTrivariateModel_SRCS src/*.C include/*.h)
if (BUILD_AS_SHARED_LIBRARY)
  ADD_LIBRARY(GoTrivariateModel SHARED ${GoTrivariateModel_SRCS})
else (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoTrivariateModel ${GoTrivariateModel_SRCS})
endif (BUILD_AS_SHARED_LIBRARY)
TARGET_LINK_LIBRARIES(GoTrivariateModel ${DEPLIBS})
SET_PROPERTY(TARGET GoTrivariateModel
  PROPERTY FOLDER "GoTrivariateModel/Libs")
SET_TARGET_PROPERTIES(GoTrivariateModel PROPERTIES SOVERSION ${GoTools_ABI_VERSION})


# Apps and tests
MACRO(ADD_APPS SUBDIR PROPERTY_FOLDER IS_TEST)
  FILE(GLOB_RECURSE GoTrivariateModel_APPS ${SUBDIR}/*.C)
  FOREACH(app ${GoTrivariateModel_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    # MESSAGE("app: " ${app})
    # MESSAGE("DEPLIBS: " ${DEPLIBS})
    TARGET_LINK_LIBRARIES(${appname} GoTrivariateModel ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${SUBDIR})
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoTrivariateModel/${PROPERTY_FOLDER}")
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
  SET(DEPLIBS ${DEPLIBS} ${Boost_LIBRARIES})
  ADD_APPS(test/unit "Unit Tests" TRUE)
  ADD_APPS(test/integration "Integration Tests" TRUE)
  ADD_APPS(test/acceptance "Acceptance Tests" TRUE)
ENDIF(GoTools_COMPILE_TESTS)


# Copy data
if (GoTools_COPY_DATA)
  FILE(COPY ${GoTrivariateModel_SOURCE_DIR}/../gotools-data/trivariatemodel/data
    DESTINATION ${GoTrivariateModel_BINARY_DIR})
  FILE(COPY ${GoTrivariateModel_SOURCE_DIR}/../gotools-data/trivariatemodel/examples/data
    DESTINATION ${GoTrivariateModel_BINARY_DIR}/examples)
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

# Install the GoTrivariateModel target with export set
install(TARGETS GoTrivariateModel
    EXPORT GoTrivariateModelTargets
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
install(EXPORT GoTrivariateModelTargets
    FILE GoTrivariateModelTargets.cmake
    #NAMESPACE GoTools::         # So consumer can do target_link_libraries(... GoTools::GoToolsCore)
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GoTrivariateModel
)

# Generate and install the Config and Version files
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/GoTrivariateModelConfigVersion.cmake"
    VERSION ${GoTools_VERSION_MAJOR}.${GoTools_VERSION_MINOR}.${GoTools_VERSION_PATCH}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/GoTrivariateModelConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/GoTrivariateModelConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GoTrivariateModel
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/GoTrivariateModelConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/GoTrivariateModelConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GoTrivariateModel
)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libGoTrivariateModel.so)
