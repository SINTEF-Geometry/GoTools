PROJECT(GoTools)


CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

# Set version info
SET(GoTools_VERSION_MAJOR 4)
SET(GoTools_VERSION_MINOR 3)
SET(GoTools_VERSION_PATCH 0)
SET(GoTools_VERSION "${GoTools_VERSION_MAJOR}.${GoTools_VERSION_MINOR}.${GoTools_VERSION_PATCH}")
SET(GoTools_ABI_VERSION "${GoTools_VERSION_MAJOR}.${GoTools_VERSION_MINOR}" CACHE INTERNAL "GoTools ABI version")

option(GoTools_COPY_DATA "Copy data?" ON)

OPTION(GoTools_ENABLE_OPENMP "Include OPENMP?" OFF)
IF(GoTools_ENABLE_OPENMP)
  FIND_PACKAGE(OpenMP REQUIRED)
ENDIF()

OPTION(GoTools_ENABLE_LOGGER "Enable logger?" OFF)
if (GoTools_ENABLE_LOGGER)
  add_definitions(-DGOTOOLS_LOG)
endif()

# Generate header with version info
#CONFIGURE_FILE(gotools-core/include/GoTools/geometry/GoTools_version.h.in
#               ${PROJECT_SOURCE_DIR}/gotools-core/include/GoTools/geometry/GoTools_version.h @ONLY)


# Set various compiler flags
ENABLE_LANGUAGE(CXX)
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wconversion -Wno-unused-but-set-variable -fPIC")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)
IF(WIN32)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP8 /DNOMINMAX /DNOCRYPT -D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS")
ENDIF(WIN32)
# Get rid of MESSAGEs in Release mode
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNVERBOSE")

# set(Boost_USE_STATIC_LIBS ON)

# Check if Boost is available
SET(Boost_ADDITIONAL_VERSIONS
  "1.74" "1.74.0"
  "1.65" "1.65.0"
  "1.58" "1.58.0"
  "1.51" "1.51.0"
  "1.49" "1.49.0"
  "1.47" "1.47.0"
  )
FIND_PACKAGE(Boost)

# Check if compiler supports c++-0x
INCLUDE(CheckCXXCompilerFlag)
IF(CMAKE_CXX_COMPILER_ID MATCHES GNU OR
   CMAKE_CXX_COMPILER_ID MATCHES Clang)
  CHECK_CXX_COMPILER_FLAG("-std=gnu++0x" HAVE_0x)
  IF(HAVE_0x)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
  ELSE(HAVE_0x)
    # C++0x is not supported - check for Boost?
    IF(Boost_FOUND)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_BOOST=1")
      SET(GoTools_COMMON_INCLUDE_DIRS
        ${GoTools_COMMON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
    ELSE(Boost_FOUND)
      MESSAGE(FATAL_ERROR "Either Boost or a compiler with c++0x support is needed")
    ENDIF(Boost_FOUND)
  ENDIF(HAVE_0x)
ENDIF() #CMAKE_CXX_COMPILER_ID MATCHES GNU or Clang)
IF(CMAKE_CXX_COMPILER_ID MATCHES Intel)
  # icpc's c++0x is lacking
  IF(Boost_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_BOOST=1")
    SET(GoTools_COMMON_INCLUDE_DIRS
      ${GoTools_COMMON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
  ELSE(Boost_FOUND)
    MESSAGE(FATAL_ERROR "Either Boost or a compiler with c++0x support is needed")
  ENDIF(Boost_FOUND)
ENDIF(CMAKE_CXX_COMPILER_ID MATCHES Intel)
IF(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
  IF(MSVC90) # Visual Studio 2008 support (c++0x is missing)
    IF(Boost_FOUND)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_BOOST=1")
      SET(GoTools_COMMON_INCLUDE_DIRS
      ${GoTools_COMMON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
    ELSE(Boost_FOUND)
      MESSAGE(FATAL_ERROR "Either Boost or a compiler with c++0x support is needed")
    ENDIF(Boost_FOUND)
  ENDIF(MSVC90)
ENDIF(CMAKE_CXX_COMPILER_ID MATCHES MSVC)


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


# Compile apps and examples
OPTION(GoTools_COMPILE_APPS
       "Compile applications?" ON)


# Compile unit test - depends on Boost
IF(Boost_FOUND)
  OPTION(GoTools_COMPILE_TESTS
    "Compile unit tests?" ON)
  IF(GoTools_COMPILE_TESTS)
    SET(Boost_USE_STATIC_LIBS ON)
    FIND_PACKAGE(Boost COMPONENTS unit_test_framework REQUIRED)
    SET(GoTools_COMMON_INCLUDE_DIRS
      ${GoTools_COMMON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
    FILE(MAKE_DIRECTORY ${GoTools_BINARY_DIR}/Testing)
    ENABLE_TESTING()
  ENDIF(GoTools_COMPILE_TESTS)
ENDIF(Boost_FOUND)


# Comment out the modules you don't want to build

ADD_SUBDIRECTORY(newmat)
ADD_SUBDIRECTORY(sisl)
ADD_SUBDIRECTORY(ttl)
SET(GoTools_COMMON_INCLUDE_DIRS
  ${newmat_SOURCE_DIR}/include
  ${sisl_SOURCE_DIR}/include
  ${ttl_SOURCE_DIR}/include
  ${GoTools_COMMON_INCLUDE_DIRS} 
  CACHE INTERNAL "GoTools common include dirs")

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
IF(GoTools_COMPILE_MODULE_qualitymodule)
ADD_SUBDIRECTORY(qualitymodule)
ENDIF(GoTools_COMPILE_MODULE_qualitymodule)

OPTION(GoTools_COMPILE_MODULE_isogeometric_model
  "Compile the GoTools module isogeometric_model?" ON)
IF(GoTools_COMPILE_MODULE_isogeometric_model)
ADD_SUBDIRECTORY(isogeometric_model)
ENDIF(GoTools_COMPILE_MODULE_isogeometric_model)

OPTION(GoTools_COMPILE_MODULE_lrsplines2D
  "Compile the GoTools module lrsplines2D?" ON)
IF(GoTools_COMPILE_MODULE_lrsplines2D)
ADD_SUBDIRECTORY(lrsplines2D)
ENDIF(GoTools_COMPILE_MODULE_lrsplines2D)

OPTION(GoTools_COMPILE_MODULE_lrsplines3D
  "Compile the GoTools module lrsplines3D?" ON)
IF(GoTools_COMPILE_MODULE_lrsplines3D)
  ADD_SUBDIRECTORY(lrsplines3D)
ENDIF(GoTools_COMPILE_MODULE_lrsplines3D)

OPTION(GoTools_COMPILE_MODULE_viewlib
  "Compile the GoTools module viewlib?" ON)
IF(GoTools_COMPILE_MODULE_viewlib)
ADD_SUBDIRECTORY(viewlib)
ENDIF(GoTools_COMPILE_MODULE_viewlib)

# CPack stuff
SET(CPACK_PACKAGE_NAME "libgotools")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "GoTools-${GoTools_VERSION}")
SET(CPACK_SOURCE_IGNORE_FILES "/CVS/;/.svn/;.swp$;.#;/#;.*~")
SET(CPACK_PACKAGE_FILE_NAME "libgotools_${CPACK_PACKAGE_VERSION}")
SET(CPACK_SOURCE_GENERATOR "ZIP;TGZ")

INCLUDE(CPack)
