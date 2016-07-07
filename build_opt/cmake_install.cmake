# Install script for directory: /home/vsk/projects/gotools_all2

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/home/vsk/projects/gotools_all2/build_opt/newmat/cmake_install.cmake")
  INCLUDE("/home/vsk/projects/gotools_all2/build_opt/sisl/cmake_install.cmake")
  INCLUDE("/home/vsk/projects/gotools_all2/build_opt/ttl/cmake_install.cmake")
  INCLUDE("/home/vsk/projects/gotools_all2/build_opt/gotools-core/cmake_install.cmake")
  INCLUDE("/home/vsk/projects/gotools_all2/build_opt/igeslib/cmake_install.cmake")
  INCLUDE("/home/vsk/projects/gotools_all2/build_opt/parametrization/cmake_install.cmake")
  INCLUDE("/home/vsk/projects/gotools_all2/build_opt/implicitization/cmake_install.cmake")
  INCLUDE("/home/vsk/projects/gotools_all2/build_opt/intersections/cmake_install.cmake")
  INCLUDE("/home/vsk/projects/gotools_all2/build_opt/trivariate/cmake_install.cmake")
  INCLUDE("/home/vsk/projects/gotools_all2/build_opt/topology/cmake_install.cmake")
  INCLUDE("/home/vsk/projects/gotools_all2/build_opt/compositemodel/cmake_install.cmake")
  INCLUDE("/home/vsk/projects/gotools_all2/build_opt/trivariatemodel/cmake_install.cmake")
  INCLUDE("/home/vsk/projects/gotools_all2/build_opt/qualitymodule/cmake_install.cmake")
  INCLUDE("/home/vsk/projects/gotools_all2/build_opt/isogeometric_model/cmake_install.cmake")
  INCLUDE("/home/vsk/projects/gotools_all2/build_opt/lrsplines2D/cmake_install.cmake")
  INCLUDE("/home/vsk/projects/gotools_all2/build_opt/viewlib/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/home/vsk/projects/gotools_all2/build_opt/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/home/vsk/projects/gotools_all2/build_opt/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
