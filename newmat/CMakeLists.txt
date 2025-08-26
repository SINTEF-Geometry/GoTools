PROJECT(newmat)


# Include directories

INCLUDE_DIRECTORIES(
  ${newmat_SOURCE_DIR}/include
  )


# Make the newmat library

FILE(GLOB_RECURSE newmat_SRCS src/*.cpp include/*.h)
if (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(newmat SHARED ${newmat_SRCS})
else (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(newmat ${newmat_SRCS})
endif (BUILD_AS_SHARED_LIBRARY)
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

# 'install' target for newmat

# Set the installation prefix based on platform
IF(WIN32)
  SET(newmat_INSTALL_PREFIX "$ENV{PROGRAMFILES}/SINTEF/newmat" CACHE PATH "Path to install newmat")
ELSE()
  SET(newmat_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Path to install newmat")
ENDIF()

# === Install Targets and Export CMake Package ===
include(CMakePackageConfigHelpers)

# Install the newmat target with export
install(TARGETS newmat
    EXPORT NewmatTargets
    ARCHIVE DESTINATION ${newmat_INSTALL_PREFIX}/lib
    LIBRARY DESTINATION ${newmat_INSTALL_PREFIX}/lib
    RUNTIME DESTINATION ${newmat_INSTALL_PREFIX}/bin
    INCLUDES DESTINATION ${newmat_INSTALL_PREFIX}/include
)

# Install header files
install(DIRECTORY include/
    DESTINATION ${newmat_INSTALL_PREFIX}/include/newmat
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
)

# Export target info to a file
install(EXPORT NewmatTargets
    FILE NewmatTargets.cmake
    #NAMESPACE Newmat::           # Uncomment if you want to use target namespaces
    DESTINATION ${newmat_INSTALL_PREFIX}/lib/cmake/newmat
)

# Generate and install version/config files
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/newmatConfigVersion.cmake"
    VERSION 10.0.0
    COMPATIBILITY SameMajorVersion
)

set(PACKAGE_INCLUDE_INSTALL_DIR "${newmat_INSTALL_PREFIX}/include")

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/newmatConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/newmatConfig.cmake"
    INSTALL_DESTINATION ${newmat_INSTALL_PREFIX}/lib/cmake/newmat
    PATH_VARS PACKAGE_INCLUDE_INSTALL_DIR
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/newmatConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/newmatConfigVersion.cmake"
    DESTINATION ${newmat_INSTALL_PREFIX}/lib/cmake/newmat
)

# CPack setup for packaging if needed
IF(NOT GoTools_ALL_COMPILE)
  SET(CPACK_SOURCE_PACKAGE_FILE_NAME "newmat-10D")
  SET(CPACK_SOURCE_IGNORE_FILES "/CVS/;/\\.svn/;\\.swp$;\\.#;/#;.*~")
  INCLUDE(CPack)
ENDIF()

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libnewmat.so)
