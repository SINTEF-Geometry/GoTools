PROJECT(GoLRspline3D)


#SET(ENABLE_OPENMP 1)

IF(GoTools_ENABLE_OPENMP)
  FIND_PACKAGE(OpenMP REQUIRED)
ENDIF(GoTools_ENABLE_OPENMP)

# Include directories

INCLUDE_DIRECTORIES(
  ${GoLRspline2D_SOURCE_DIR}/include
  ${GoLRspline3D_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${GoTrivariate_SOURCE_DIR}/include
  ${GoTools_COMMON_INCLUDE_DIRS}
  )


# Linked in libraries

SET(DEPLIBS
  GoToolsCore
  sisl
  GoTrivariate
  GoLRspline2D
  )

#message(STATUS DEBLIPS:${DEPLIBS})

# Make the GoLRspline3D library

FILE(GLOB_RECURSE GoLRspline3D_SRCS src/*.C include/*.h)
if (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoLRspline3D SHARED ${GoLRspline3D_SRCS})
else (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoLRspline3D ${GoLRspline3D_SRCS})
endif (BUILD_AS_SHARED_LIBRARY)
TARGET_LINK_LIBRARIES(GoLRspline3D ${DEPLIBS})
SET_PROPERTY(TARGET GoLRspline3D
  PROPERTY FOLDER "GoLRspline3D/Libs")
SET_TARGET_PROPERTIES(GoLRspline3D PROPERTIES SOVERSION ${GoTools_ABI_VERSION})
IF(GoTools_ENABLE_OPENMP)
  SET_TARGET_PROPERTIES(GoLRspline3D PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}") 
  SET_TARGET_PROPERTIES(GoLRspline3D PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
ENDIF(GoTools_ENABLE_OPENMP)

#message(STATUS ${GoLRspline3D_SRCS})


# # Make the LRSplinePlotUtils library

# FILE(GLOB_RECURSE LRspline3D_SRCS src/*.C include/*.h)
# ADD_LIBRARY(GoLRspline3D ${GoLRspline3D_SRCS})
# TARGET_LINK_LIBRARIES(GoLRspline3D ${DEPLIBS})
# SET_PROPERTY(TARGET GoLRspline3D
#   PROPERTY FOLDER "GoLRspline3D/Libs")


# Apps, examples, tests, ...?

# Apps and tests
MACRO(ADD_APPS SUBDIR PROPERTY_FOLDER IS_TEST)
  FILE(GLOB_RECURSE GoLRspline3D_APPS ${SUBDIR}/*.C)
  FOREACH(app ${GoLRspline3D_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoLRspline3D ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${SUBDIR})
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoLRspline3D/${PROPERTY_FOLDER}")
    IF(${IS_TEST})
      ADD_TEST(${appname} ${SUBDIR}/${appname}
	--log_format=XML --log_level=all --log_sink=../Testing/${appname}.xml)
      SET_TESTS_PROPERTIES( ${appname} PROPERTIES LABELS "${SUBDIR}" )
    ENDIF(${IS_TEST})
  ENDFOREACH(app)
ENDMACRO(ADD_APPS)

IF(GoTools_COMPILE_APPS)
  FILE(GLOB_RECURSE GoLRspline3D_APPS app/*.C)
  FOREACH(app ${GoLRspline3D_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
#    message(STATUS ${appname})
#    message(STATUS ${app})
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoLRspline3D ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
    IF(GoTools_ENABLE_OPENMP)
      SET_TARGET_PROPERTIES(${appname} PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}") 
      SET_TARGET_PROPERTIES(${appname} PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
    ENDIF(GoTools_ENABLE_OPENMP)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoLRspline3D/Apps")
  ENDFOREACH(app)

  FILE(GLOB_RECURSE GoLRspline3D_EXAMPLES examples/*.C)
  FOREACH(app ${GoLRspline3D_EXAMPLES})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoLRspline3D ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY examples)
    IF(GoTools_ENABLE_OPENMP)
      SET_TARGET_PROPERTIES(${appname} PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}") 
      SET_TARGET_PROPERTIES(${appname} PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
    ENDIF(GoTools_ENABLE_OPENMP)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoLRspline3D/Examples")
  ENDFOREACH(app)
ENDIF(GoTools_COMPILE_APPS)

IF(GoTools_COMPILE_TESTS)
  SET(DEPLIBS ${DEPLIBS} ${Boost_LIBRARIES})
  ADD_APPS(test/unit "Unit Tests" TRUE)
  ADD_APPS(test/integration "Integration Tests" TRUE)
  ADD_APPS(test/acceptance "Acceptance Tests" TRUE)
ENDIF(GoTools_COMPILE_TESTS)

# Copy/symlink data dir.
IF(GoTools_COPY_DATA)
  if (WIN32)
    FILE(COPY ${GoLRspline3D_SOURCE_DIR}/data
      DESTINATION ${GoLRspline3D_BINARY_DIR})
  else()
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data
      ${CMAKE_CURRENT_BINARY_DIR}/data)
  endif()
ENDIF()

# Create a tmp directory.
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tmp)

# 'install' target

IF(WIN32)
  # Windows
  # lib
  INSTALL(TARGETS GoLRspline3D DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  # include
  INSTALL(DIRECTORY include/GoTools/lrsplines3D
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS GoLRspline3D DESTINATION lib COMPONENT lrsplines3D)
  # include
  INSTALL(DIRECTORY include/GoTools/lrsplines3D
    DESTINATION include/GoTools
    COMPONENT lrsplines3D-dev
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libGoLRspline3D.so)
