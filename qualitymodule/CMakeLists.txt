PROJECT(GoQualityModule)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoQualityModule_SOURCE_DIR}/include
  ${GoCompositeModel_SOURCE_DIR}/include
  ${parametrization_SOURCE_DIR}/include
  ${GoTopology_SOURCE_DIR}/include
  ${GoIntersections_SOURCE_DIR}/include
  ${GoImplicitization_SOURCE_DIR}/include
  ${GoIgeslib_SOURCE_DIR}/include
  ${GoTrivariate_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${sisl_SOURCE_DIR}/include
  ${newmat_SOURCE_DIR}/include
  )


# Linked in libraries

SET(DEPLIBS
  GoCompositeModel
  parametrization
  GoTopology
  GoIntersections
  GoImplicitization
  GoIgeslib
  GoToolsCore
  sisl
  newmat
  )


# Make the GoQualityModule library

FILE(GLOB_RECURSE GoQualityModule_SRCS src/*.C include/*.h)
ADD_LIBRARY(GoQualityModule ${GoQualityModule_SRCS})
#TARGET_LINK_LIBRARIES(GoQualityModule ${DEPLIBS})


# Apps, examples, tests, ...?
IF(GoTools_COMPILE_APPS)
  FILE(GLOB GoQualityModule_APPS app/*.C)
  FOREACH(app ${GoQualityModule_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoQualityModule ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
  ENDFOREACH(app)
ENDIF(GoTools_COMPILE_APPS)

# 'install' target

IF(WIN32)
  # Windows
  # lib
  INSTALL(TARGETS GoQualityModule DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  # include
  INSTALL(DIRECTORY include/GoTools/qualitymodule
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS GoQualityModule DESTINATION lib/GoTools)
  # include
  INSTALL(DIRECTORY include/GoTools/qualitymodule
    DESTINATION include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)
