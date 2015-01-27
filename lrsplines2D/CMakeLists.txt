PROJECT(GoLRspline2D)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoLRspline2D_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${GoTools_COMMON_INCLUDE_DIRS}
  )


# Linked in libraries

SET(DEPLIBS
  GoToolsCore
  sisl
  )


# Make the GoLRspline2D library

FILE(GLOB_RECURSE GoLRspline2D_SRCS src/*.C include/*.h)
ADD_LIBRARY(GoLRspline2D ${GoLRspline2D_SRCS})
TARGET_LINK_LIBRARIES(GoLRspline2D ${DEPLIBS})
SET_PROPERTY(TARGET GoLRspline2D
  PROPERTY FOLDER "GoLRspline2D/Libs")
SET_TARGET_PROPERTIES(GoLRspline2D PROPERTIES SOVERSION ${GoTools_ABI_VERSION})

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
    ENDIF(${IS_TEST})
  ENDFOREACH(app)
ENDMACRO(ADD_APPS)

IF(GoTools_COMPILE_APPS)
  FILE(GLOB_RECURSE GoLRspline2D_APPS app/*.C)
  FOREACH(app ${GoLRspline2D_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
#    message(STATUS ${appname})
#    message(STATUS ${app})
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoLRspline2D ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
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
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoLRspline2D/Examples")
  ENDFOREACH(app)
ENDIF(GoTools_COMPILE_APPS)


IF(GoTools_COMPILE_TESTS)
  SET(DEPLIBS ${DEPLIBS} ${Boost_LIBRARIES})
  ADD_APPS(unit_tests "Unit Tests" TRUE)
  ADD_APPS(regression_tests "Regression Tests" TRUE)
ENDIF(GoTools_COMPILE_TESTS)

# Symlink data dir.
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data
  ${CMAKE_CURRENT_BINARY_DIR}/data)

# Create a tmp directory.
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tmp)

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
