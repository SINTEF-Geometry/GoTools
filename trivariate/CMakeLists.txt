PROJECT(GoTrivariate)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoTrivariate_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  #  ${sisl_SOURCE_DIR}/include
  )

IF(USE_BOOST)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_BOOST=1")
  FIND_PACKAGE(Boost REQUIRED)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ELSE(USE_BOOST)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
ENDIF(USE_BOOST)

# Linked in libraries

SET(DEPLIBS
  GoToolsCore
  #  sisl
  )

SET(DEPLIBS "${DEPLIBS} -L/usr/local/lib/GoTools")

# Make the GoTrivariate library

FILE(GLOB_RECURSE GoTrivariate_SRCS src/*.C include/*.h)
ADD_LIBRARY(GoTrivariate ${GoTrivariate_SRCS})


# Apps, examples, tests, ...?

#FILE(GLOB GoTrivariate_APPS app/*.C)
#FOREACH(app ${GoTrivariate_APPS})
#  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
#  ADD_EXECUTABLE(${appname} ${app})
#  TARGET_LINK_LIBRARIES(${appname} GoTrivariate ${DEPLIBS})
#ENDFOREACH(app)


# 'install' target

IF(WIN32)
  # Windows
  # lib
  INSTALL(TARGETS GoTrivariate DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  # include
  INSTALL(DIRECTORY include/GoTools/trivariate
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS GoTrivariate DESTINATION lib/GoTools)
  # include
  INSTALL(DIRECTORY include/GoTools/trivariate
    DESTINATION include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)
