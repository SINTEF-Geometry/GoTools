PROJECT(GoTrivariateModel)

# Find modules

#FIND_PACKAGE(PugiXML REQUIRED)

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
  ${PUGIXML_LIBRARY}
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
    ENDIF(${IS_TEST})
  ENDFOREACH(app)
ENDMACRO(ADD_APPS)

IF(GoTools_COMPILE_APPS)
  ADD_APPS(app "Apps" FALSE)
  ADD_APPS(examples "Examples" FALSE)
ENDIF(GoTools_COMPILE_APPS)

IF(GoTools_COMPILE_TESTS)
  SET(DEPLIBS ${DEPLIBS} ${Boost_LIBRARIES})
  ADD_APPS(unit_tests "Unit Tests" TRUE)
  ADD_APPS(regression_tests "Regression Tests" TRUE)
ENDIF(GoTools_COMPILE_TESTS)


# Copy data
FILE(COPY ${GoTrivariateModel_SOURCE_DIR}/data
  DESTINATION ${GoTrivariateModel_BINARY_DIR})
FILE(COPY ${GoTrivariateModel_SOURCE_DIR}/examples/data
  DESTINATION ${GoTrivariateModel_BINARY_DIR}/examples)


# 'install' target

IF(WIN32)
  # Windows
  # lib
  INSTALL(TARGETS GoTrivariateModel DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  #include
  INSTALL(DIRECTORY include/GoTools/trivariatemodel
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  #lib
  INSTALL(TARGETS GoTrivariateModel DESTINATION lib COMPONENT trivariatemodel)
  # include
  INSTALL(DIRECTORY include/GoTools/trivariatemodel
    DESTINATION include/GoTools
    COMPONENT trivariatemodel-dev
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libGoTrivariateModel.so)
