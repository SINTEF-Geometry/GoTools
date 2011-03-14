PROJECT(newmat)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


# Include directories

INCLUDE_DIRECTORIES(
  ${newmat_SOURCE_DIR}/include
  )


# Make the newmat library

FILE(GLOB newmat_SRCS src/*.cpp)
ADD_LIBRARY(newmat ${newmat_SRCS})


# Apps, examples, tests, ...?

FILE(GLOB tmt_SRCS app/tmt.cpp app/tmt?.cpp)
ADD_EXECUTABLE(tmt ${tmt_SRCS})
TARGET_LINK_LIBRARIES(tmt newmat)

ADD_EXECUTABLE(example app/example.cpp)
TARGET_LINK_LIBRARIES(example newmat)

ADD_EXECUTABLE(garch app/garch.cpp)
TARGET_LINK_LIBRARIES(garch newmat)
ADD_CUSTOM_COMMAND(
  TARGET garch
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${newmat_SOURCE_DIR}/data/garch.dat ${newmat_BINARY_DIR}
  )

ADD_EXECUTABLE(nl_ex app/nl_ex.cpp)
TARGET_LINK_LIBRARIES(nl_ex newmat)

ADD_EXECUTABLE(sl_ex app/sl_ex.cpp)
TARGET_LINK_LIBRARIES(sl_ex newmat)

ADD_EXECUTABLE(test_exc app/test_exc.cpp)
TARGET_LINK_LIBRARIES(test_exc newmat)


# Install newmat
INSTALL(TARGETS newmat DESTINATION lib)
IF(WIN32)
  # Windows
  INSTALL(DIRECTORY include/
    DESTINATION include/newmat
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  INSTALL(DIRECTORY include/
    DESTINATION include/newmat
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)