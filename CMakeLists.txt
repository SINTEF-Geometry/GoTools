PROJECT(GoTools)


CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(GoTools_VERSION_MAJOR 4)
SET(GoTools_VERSION_MINOR 0)
SET(GoTools_VERSION_PATCH 1)
SET(GoTools_VERSION "${GoTools_VERSION_MAJOR}.${GoTools_VERSION_MINOR}.${GoTools_VERSION_PATCH}")

# Generate header with version info
CONFIGURE_FILE(gotools-core/include/GoTools/geometry/GoTools_version.h.in
               ${PROJECT_SOURCE_DIR}/gotools-core/include/GoTools/geometry/GoTools_version.h @ONLY)
CONFIGURE_FILE(gpl/fixupdebs.sh.in
               ${PROJECT_SOURCE_DIR}/gpl/fixupdebs.sh @ONLY)
CONFIGURE_FILE(gpl/generate_gpl_headers.in
               ${PROJECT_SOURCE_DIR}/gpl/generate_gpl_headers @ONLY)

ENABLE_LANGUAGE(CXX)

# Set various compiler flags
IF(NOT WIN32)
  # Set warning level to '-Wall' on Linux
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wconversion -Wno-unused-but-set-variable -fPIC")
ELSE(NOT WIN32)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP8 /DNOMINMAX /DNOCRYPT
  -D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS")
ENDIF(NOT WIN32)

# Get rid of MESSAGEs in Release mode
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNVERBOSE")

OPTION(GoTools_USE_BOOST
       "Use BOOST instead of c++0x?" OFF)

OPTION(GoTools_COMPILE_APPS
       "Compile applications?" ON)
IF(GoTools_USE_BOOST)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_BOOST=1")
  FIND_PACKAGE(Boost REQUIRED)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ELSE(GoTools_USE_BOOST)
  INCLUDE(CheckCXXCompilerFlag)
  IF(CMAKE_CXX_COMPILER_ID MATCHES GNU)
  # check if compiler supports c++-0x
    CHECK_CXX_COMPILER_FLAG("-std=gnu++0x" HAVE_0x)
    IF(HAVE_0x)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
    ELSE(HAVE_0x)
      MESSAGE(FATAL_ERROR "A compiler with c++-0x support is needed")
    ENDIF(HAVE_0x)
  ENDIF(CMAKE_CXX_COMPILER_ID MATCHES GNU)
ENDIF(GoTools_USE_BOOST)


# Set install prefix on Windows
IF(WIN32)
  SET(CMAKE_INSTALL_PREFIX CACHE INTERNAL "")
  SET(GoTools_INSTALL_PREFIX "$ENV{PROGRAMFILES}/SINTEF/GoTools"
    CACHE PATH "Path to install GoTools")
ENDIF(WIN32)


# Organize the project in folders (VS only?)
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)

# For newmat/sisl/ttl to know they are compiled as modules

SET(GoTools_ALL_COMPILE 1)


# Comment out the modules you don't want to build

OPTION(GoTools_COMPILE_GPL
       "Compile the GPL modules only?" OFF)

ADD_SUBDIRECTORY(newmat)
ADD_SUBDIRECTORY(sisl)
ADD_SUBDIRECTORY(ttl)

ADD_SUBDIRECTORY(gotools-core)

OPTION(GoTools_COMPILE_MODULE_igeslib 
  "Compile the GoTools module igeslib?" ON)
IF(GoTools_COMPILE_MODULE_igeslib)
ADD_SUBDIRECTORY(igeslib)
ENDIF(GoTools_COMPILE_MODULE_igeslib)

OPTION(GoTools_COMPILE_MODULE_parametrization 
  "Compile the GoTools module parametrization?" ON)
IF(GoTools_COMPILE_MODULE_parametrization)
ADD_SUBDIRECTORY(parametrization)
ENDIF(GoTools_COMPILE_MODULE_parametrization)

OPTION(GoTools_COMPILE_MODULE_implicitization
  "Compile the GoTools module implicitization?" ON)
IF(GoTools_COMPILE_MODULE_implicitization)
ADD_SUBDIRECTORY(implicitization)
ENDIF(GoTools_COMPILE_MODULE_implicitization)

OPTION(GoTools_COMPILE_MODULE_intersections
  "Compile the GoTools module intersections?" ON)
IF(GoTools_COMPILE_MODULE_intersections)
ADD_SUBDIRECTORY(intersections)
ENDIF(GoTools_COMPILE_MODULE_intersections)

OPTION(GoTools_COMPILE_MODULE_trivariate
  "Compile the GoTools module trivariate?" ON)
IF(GoTools_COMPILE_MODULE_trivariate)
ADD_SUBDIRECTORY(trivariate)
ENDIF(GoTools_COMPILE_MODULE_trivariate)

OPTION(GoTools_COMPILE_MODULE_viewlib
  "Compile the GoTools module viewlib?" ON)
IF(GoTools_COMPILE_MODULE_viewlib)
ADD_SUBDIRECTORY(viewlib)
ENDIF(GoTools_COMPILE_MODULE_viewlib)

OPTION(GoTools_COMPILE_MODULE_topology
  "Compile the GoTools module topology?" ON)
IF(GoTools_COMPILE_MODULE_topology)
ADD_SUBDIRECTORY(topology)
ENDIF(GoTools_COMPILE_MODULE_topology)

OPTION(GoTools_COMPILE_MODULE_compositemodel
  "Compile the GoTools module compositemodel?" ON)
IF(GoTools_COMPILE_MODULE_compositemodel)
ADD_SUBDIRECTORY(compositemodel)
ENDIF(GoTools_COMPILE_MODULE_compositemodel)

OPTION(GoTools_COMPILE_MODULE_trivariatemodel
  "Compile the GoTools module trivariatemodel?" ON)
IF(GoTools_COMPILE_MODULE_trivariatemodel)
ADD_SUBDIRECTORY(trivariatemodel)
ENDIF(GoTools_COMPILE_MODULE_trivariatemodel)

OPTION(GoTools_COMPILE_MODULE_qualitymodule
  "Compile the GoTools module qualitymodule?" ON)
IF(GoTools_COMPILE_MODULE_qualitymodule AND NOT GoTools_COMPILE_GPL)
ADD_SUBDIRECTORY(qualitymodule)
ENDIF(GoTools_COMPILE_MODULE_qualitymodule AND NOT GoTools_COMPILE_GPL)

OPTION(GoTools_COMPILE_MODULE_model_toolbox
  "Compile the GoTools modulemodel_toolbox ?" ON)
IF(GoTools_COMPILE_MODULE_model_toolbox AND NOT GoTools_COMPILE_GPL)
ADD_SUBDIRECTORY(model_toolbox)
ENDIF(GoTools_COMPILE_MODULE_model_toolbox AND NOT GoTools_COMPILE_GPL)

OPTION(GoTools_COMPILE_MODULE_isogeometric_model
  "Compile the GoTools module isogeometric_model?" ON)
IF(GoTools_COMPILE_MODULE_isogeometric_model AND NOT GoTools_COMPILE_GPL)
ADD_SUBDIRECTORY(isogeometric_model)
ENDIF(GoTools_COMPILE_MODULE_isogeometric_model AND NOT GoTools_COMPILE_GPL)

OPTION(GoTools_COMPILE_MODULE_lr_splines
  "Compile the GoTools module lr_splines? Depends on Boost." OFF)
IF(GoTools_COMPILE_MODULE_lr_splines AND NOT GoTools_COMPILE_GPL)
ADD_SUBDIRECTORY(lr_splines)
ENDIF(GoTools_COMPILE_MODULE_lr_splines AND NOT GoTools_COMPILE_GPL)

OPTION(GoTools_COMPILE_MODULE_geometric_modeler
  "Compile the GoTools module geometric_modeler? " OFF)
IF(GoTools_COMPILE_MODULE_geometric_modeler AND NOT GoTools_COMPILE_GPL)
ADD_SUBDIRECTORY(geometric_modeler)
ENDIF(GoTools_COMPILE_MODULE_geometric_modeler AND NOT GoTools_COMPILE_GPL)

IF(WIN32)
  OPTION(GoTools_COMPILE_MODULE_step_reader
    "Compile the GoTools module step_reader? Depends on EPM." OFF)
  IF(GoTools_COMPILE_MODULE_step_reader)
    ADD_SUBDIRECTORY(step_reader)
  ENDIF(GoTools_COMPILE_MODULE_step_reader)

  OPTION(GoTools_COMPILE_MODULE_step_reader_cpp
    "Compile the GoTools module step_reader_cpp? Depends on EPM, C++ version."
    OFF)
  IF(GoTools_COMPILE_MODULE_step_reader_cpp)
    ADD_SUBDIRECTORY(step_reader_cpp)
  ENDIF(GoTools_COMPILE_MODULE_step_reader_cpp)
ENDIF(WIN32)



# CPack stuff
SET(CPACK_PACKAGE_NAME "libgotools")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Arne Morten Kvarving <arne.morten.kvarving@sintef.no>")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The GoTools geometry/spline library")
SET(CPACK_PACKAGE_VERSION "${GoTools_VERSION}")
SET(CPACK_PACKAGE_VERSION_MAJOR ${GoTools_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${GoTools_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${GoTools_VERSION_PATCH})
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "GoTools-${GoTools_VERSION}")
SET(CPACK_SOURCE_IGNORE_FILES "/CVS/;/.svn/;.swp$;.#;/#;.*~")
SET(CPACK_GENERATOR "DEB")
SET(CPACK_PACKAGE_FILE_NAME "libgotools_${CPACK_PACKAGE_VERSION}")
SET(CPACK_DEB_COMPONENT_INSTALL ON)
SET(CPACK_SOURCE_GENERATOR "ZIP;TGZ")
IF (NOT WIN32)
  EXECUTE_PROCESS(COMMAND lsb_release "-sc" OUTPUT_VARIABLE CODENAME)
  STRING(REGEX REPLACE "\n" "" SYSTEM_CODENAME ${CODENAME})
  EXECUTE_PROCESS(COMMAND dpkg --print-architecture OUTPUT_VARIABLE ARCHITECTURE)
  STRING(REGEX REPLACE "\n" "" CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${ARCHITECTURE})
  SET(CPACK_PACKAGE_FILE_NAME "libgotools_${CPACK_PACKAGE_VERSION}_${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}-${SYSTEM_CODENAME}")
ENDIF(NOT WIN32)

INCLUDE(CPack)

# For generating ubuntu packages                                
ADD_CUSTOM_TARGET(ubuntu 
                  COMMAND make package
                  COMMAND ${PROJECT_SOURCE_DIR}/gpl/fixupdebs.sh
                  COMMAND echo "All done. Packages are in UbuntuDebs/"
                  COMMENT "Generating and fixing up Ubuntu packages" VERBATIM)
