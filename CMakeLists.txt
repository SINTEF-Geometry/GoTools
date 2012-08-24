PROJECT(newmat)


# Include directories

INCLUDE_DIRECTORIES(
  ${newmat_SOURCE_DIR}/include
  )


# Make the newmat library

FILE(GLOB_RECURSE newmat_SRCS src/*.cpp include/*.h)
ADD_LIBRARY(newmat ${newmat_SRCS})
SET_PROPERTY(TARGET newmat
  PROPERTY FOLDER "newmat/Libs")


# Apps, examples, tests, ...?

FILE(GLOB tmt_SRCS app/tmt.cpp app/tmt?.cpp)
ADD_EXECUTABLE(tmt ${tmt_SRCS})
TARGET_LINK_LIBRARIES(tmt newmat)
SET_PROPERTY(TARGET tmt
  PROPERTY FOLDER "newmat/Apps")

ADD_EXECUTABLE(example app/example.cpp)
TARGET_LINK_LIBRARIES(example newmat)
SET_PROPERTY(TARGET example
  PROPERTY FOLDER "newmat/Apps")

ADD_EXECUTABLE(garch app/garch.cpp)
TARGET_LINK_LIBRARIES(garch newmat)
SET_PROPERTY(TARGET garch
  PROPERTY FOLDER "newmat/Apps")
ADD_CUSTOM_COMMAND(
  TARGET garch
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${newmat_SOURCE_DIR}/data/garch.dat ${newmat_BINARY_DIR}
  )

ADD_EXECUTABLE(nl_ex app/nl_ex.cpp)
TARGET_LINK_LIBRARIES(nl_ex newmat)
SET_PROPERTY(TARGET nl_ex
  PROPERTY FOLDER "newmat/Apps")

ADD_EXECUTABLE(sl_ex app/sl_ex.cpp)
TARGET_LINK_LIBRARIES(sl_ex newmat)
SET_PROPERTY(TARGET sl_ex
  PROPERTY FOLDER "newmat/Apps")

ADD_EXECUTABLE(test_exc app/test_exc.cpp)
TARGET_LINK_LIBRARIES(test_exc newmat)
SET_PROPERTY(TARGET test_exc
  PROPERTY FOLDER "newmat/Apps")


# Install newmat

IF(WIN32)
  # Windows
  SET(CMAKE_INSTALL_PREFIX CACHE INTERNAL "")
  SET(newmat_INSTALL_PREFIX "$ENV{PROGRAMFILES}/SINTEF/newmat"
    CACHE PATH "Path to install Newmat")
  INSTALL(TARGETS newmat DESTINATION ${newmat_INSTALL_PREFIX}/lib)
  INSTALL(DIRECTORY include/
    DESTINATION ${newmat_INSTALL_PREFIX}/include/newmat
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS newmat DESTINATION lib COMPONENT newmat)
  # include
  INSTALL(DIRECTORY include/
    COMPONENT newmat-dev
    DESTINATION include/newmat
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)

# CPack stuff
IF(NOT GoTools_ALL_COMPILE)
  SET(CPACK_SOURCE_PACKAGE_FILE_NAME "newmat-4.5.0")
  SET(CPACK_SOURCE_IGNORE_FILES "/CVS/;/\\.svn/;\\.swp$;\\.#;/#;.*~")
  INCLUDE(CPack)
ENDIF(NOT GoTools_ALL_COMPILE)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libnewmat.so)
