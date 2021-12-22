PROJECT(GoLRspline2D)

#SET(ENABLE_OPENMP 1)

IF(GoTools_ENABLE_OPENMP)
  FIND_PACKAGE(OpenMP REQUIRED)
ENDIF(GoTools_ENABLE_OPENMP)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoLRspline2D_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${GoTools_COMMON_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
)


# Linked in libraries

SET(DEPLIBS
  GoToolsCore
  sisl
  )

# Make the GoLRspline2D library

FILE(GLOB_RECURSE GoLRspline2D_SRCS src/*.C include/*.h)
if (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoLRspline2D SHARED ${GoLRspline2D_SRCS})
else (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoLRspline2D ${GoLRspline2D_SRCS})
endif (BUILD_AS_SHARED_LIBRARY)
TARGET_LINK_LIBRARIES(GoLRspline2D ${DEPLIBS})
SET_PROPERTY(TARGET GoLRspline2D
  PROPERTY FOLDER "GoLRspline2D/Libs")
SET_TARGET_PROPERTIES(GoLRspline2D PROPERTIES SOVERSION ${GoTools_ABI_VERSION})
IF(GoTools_ENABLE_OPENMP)
  SET_TARGET_PROPERTIES(GoLRspline2D PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}") 
  SET_TARGET_PROPERTIES(GoLRspline2D PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
ENDIF(GoTools_ENABLE_OPENMP)


#message(STATUS ${GoLRspline2D_SRCS})


# # Make the LRSplinePlotUtils library

# FILE(GLOB_RECURSE LRspline2D_SRCS src/*.C include/*.h)
# ADD_LIBRARY(GoLRspline2D ${GoLRspline2D_SRCS})
# TARGET_LINK_LIBRARIES(GoLRspline2D ${DEPLIBS})
# SET_PROPERTY(TARGET GoLRspline2D
#   PROPERTY FOLDER "GoLRspline2D/Libs")


# Apps, examples, tests, ...?

# Apps and tests
MACRO(ADD_APPS SUBDIR PROPERTY_FOLDER IS_TEST)
  FILE(GLOB_RECURSE GoLRspline2D_APPS ${SUBDIR}/*.C)
  FOREACH(app ${GoLRspline2D_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoLRspline2D ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${SUBDIR})
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoLRspline2D/${PROPERTY_FOLDER}")
    IF(${IS_TEST})
      ADD_TEST(${appname} ${SUBDIR}/${appname}
		--log_format=XML --log_level=all --log_sink=../Testing/${appname}.xml)
      SET_TESTS_PROPERTIES( ${appname} PROPERTIES LABELS "${SUBDIR}" )
    ENDIF(${IS_TEST})
  ENDFOREACH(app)
ENDMACRO(ADD_APPS)

IF(GoTools_COMPILE_APPS)
  FILE(GLOB_RECURSE GoLRspline2D_APPS app/*.C)
  FOREACH(app ${GoLRspline2D_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
   # message(STATUS ${appname})
   # message(STATUS ${app})
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoLRspline2D ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
    IF(GoTools_ENABLE_OPENMP)
      SET_TARGET_PROPERTIES(${appname} PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}") 
      SET_TARGET_PROPERTIES(${appname} PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
    ENDIF(GoTools_ENABLE_OPENMP)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoLRspline2D/Apps")
  ENDFOREACH(app)

  FILE(GLOB_RECURSE GoLRspline2D_EXAMPLES examples/*.C)
  FOREACH(app ${GoLRspline2D_EXAMPLES})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoLRspline2D ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY examples)
    IF(GoTools_ENABLE_OPENMP)
      SET_TARGET_PROPERTIES(${appname} PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}") 
      SET_TARGET_PROPERTIES(${appname} PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
    ENDIF(GoTools_ENABLE_OPENMP)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoLRspline2D/Examples")
  ENDFOREACH(app)
ENDIF(GoTools_COMPILE_APPS)


IF(GoTools_COMPILE_TESTS)
  SET(DEPLIBS ${DEPLIBS} ${Boost_LIBRARIES})
  ADD_APPS(test/unit "Unit Tests" TRUE)
  ADD_APPS(test/integration "Integration Tests" TRUE)
  ADD_APPS(test/acceptance "Acceptance Tests" TRUE)
ENDIF(GoTools_COMPILE_TESTS)

# Create a tmp directory.
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tmp)

# Copy data
if (GoTools_COPY_DATA)
  FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../gotools-data/lrsplines2D/data
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

# 'install' target

IF(WIN32)
  # Windows
  # lib
  INSTALL(TARGETS GoLRspline2D DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  # include
  INSTALL(DIRECTORY include/GoTools/lrsplines2D
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS GoLRspline2D DESTINATION lib COMPONENT lrsplines2D)
  # include
  INSTALL(DIRECTORY include/GoTools/lrsplines2D
    DESTINATION include/GoTools
    COMPONENT lrsplines2D-dev
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libGoLRspline2D.so)
