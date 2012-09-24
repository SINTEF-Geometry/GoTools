PROJECT(GoToolsCore)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoToolsCore_SOURCE_DIR}/include
  )


# Linked in libraries


# Make the gotools-core library

FILE(GLOB_RECURSE GoToolsCore_SRCS src/*.cpp src/*.C include/*.h)
ADD_LIBRARY(GoToolsCore ${GoToolsCore_SRCS})
SET_PROPERTY(TARGET GoToolsCore
  PROPERTY FOLDER "GoToolsCore/Libs")


# Apps and examples
IF(GoTools_COMPILE_APPS)
  FILE(GLOB_RECURSE gotools_APPS app/*.C)
  FOREACH(app ${gotools_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoToolsCore)
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoToolsCore/Apps")
  ENDFOREACH(app)

  FILE(GLOB_RECURSE gotools_EXAMPLES examples/*.C)
  FOREACH(app ${gotools_EXAMPLES})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoToolsCore)
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY examples)
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoToolsCore/Examples")
  ENDFOREACH(app)
ENDIF(GoTools_COMPILE_APPS)


# Unit tests
IF(GoTools_COMPILE_TESTS)
#  SET(Boost_USE_STATIC_LIBS ON)
#  FIND_PACKAGE(Boost REQUIRED COMPONENTS unit_test_framework)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  FILE(GLOB_RECURSE gotools_TESTS test/*.C)
  FOREACH(test ${gotools_TESTS})
    GET_FILENAME_COMPONENT(testname ${test} NAME_WE)
    ADD_EXECUTABLE(${testname} ${test})
    TARGET_LINK_LIBRARIES(${testname} GoToolsCore ${Boost_LIBRARIES})
#    SET_TARGET_PROPERTIES(${testname}
#      PROPERTIES RUNTIME_OUTPUT_DIRECTORY test)
    SET_PROPERTY(TARGET ${testname}
      PROPERTY FOLDER "GoToolsCore/Tests")
    ADD_TEST(${testname} ${EXECUTABLE_OUTPUT_PATH}/${testname})
  ENDFOREACH(test)
ENDIF(GoTools_COMPILE_TESTS)


# Copy data
ADD_CUSTOM_COMMAND(
  TARGET GoToolsCore
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy_directory ${GoToolsCore_SOURCE_DIR}/data
  ${GoToolsCore_BINARY_DIR}/data
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy_directory ${GoToolsCore_SOURCE_DIR}/examples/data
  ${GoToolsCore_BINARY_DIR}/examples/data
  )


# 'install' target
IF(WIN32)
  # Windows
  # lib
  INSTALL(TARGETS GoToolsCore DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  # include
  INSTALL(DIRECTORY include/
    DESTINATION ${GoTools_INSTALL_PREFIX}/include
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS GoToolsCore DESTINATION lib COMPONENT core)
  # include
  INSTALL(DIRECTORY include/
    DESTINATION include
    COMPONENT core-dev
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libGoToolsCore.so)
